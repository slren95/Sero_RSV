library(shiny)
library(bslib)
library(tidyverse)
library(reactable)
library(sparkline)
library(plotly)
library(DT)

# 1. 读取数据
df <- readRDS('rda/data.1.rds')

# 统一预处理：将缺失值转化为可读的文本，确保筛选器能正常解析
df <- df %>%
  mutate(
    red_mark = ifelse(is.na(red_mark) | red_mark == "", "Normal", as.character(red_mark)),
    age = as.character(age),
    group = as.character(group)
  )

metrics_cols <- c(
  "ig_g_net_mfi_new_beads", "ig_a_net_mfi_new_beads",
  "ig_a_normalized_net_mfi", "ig_g_normalized_net_mfi",
  "site_iii_mpe8_inhibition", "site_ii_motavizumab_inhibition",
  "site_v_suptavumab_inhibition", "site_ph_nirsevimab_inhibition",
  "site_iv_clesrovimab_inhibition"
)

# 2. 数据预处理：提前为面板 1 做好聚合转换（保持全员静态或全数据状态）
df_all_summary <- df %>%
  group_by(subject_ID, age, group) %>%
  summarise(
    across(all_of(metrics_cols), ~ list(.x), .names = "{col}_trend"),
    .groups = "drop"
  )

# 3. UI 界面
ui <- page_navbar(
  theme = bs_theme(version = 5, bootswatch = "cosmo"),
  title = "受试者抗体科研分析平台",
  id = "tabs",

  # 面板 1：趋势纵览（Reactable + Sparkline）
  nav_panel(
    title = "趋势纵览 (Sparklines)",
    value = "summary",
    card(
      card_header("受试者抗体趋势纵览 (Visit 1-5 全员合并)"),
      reactableOutput("summary_table"),
      full_screen = TRUE
    )
  ),

  # 面板 2：明细与交互曲线（联合筛选 + DT多选 + 原生Plotly）
  nav_panel(
    title = "明细与交互曲线 (Plotly)",
    value = "detail",

    page_sidebar(
      fillable = FALSE,
      # 升级侧边栏：支持 Age, Group, Red Mark 三重联合筛选
      sidebar = sidebar(
        title = "组合筛选面板",

        checkboxGroupInput(
          "selected_ages",
          "1. 筛选年龄 (Age):",
          choices = sort(unique(df$age)),
          selected = unique(df$age)
        ),
        hr(),
        checkboxGroupInput(
          "selected_groups",
          "2. 筛选年龄组 (Group):",
          choices = sort(unique(df$group)),
          selected = unique(df$group)
        ),
        hr(),
        checkboxGroupInput(
          "selected_marks",
          "3. 筛选红标状态 (Red Mark):",
          choices = sort(unique(df$red_mark)),
          selected = unique(df$red_mark)
        )
      ),

      fluidRow(
        # 左侧 4格：使用 DT 展示包含红标的唯一个体清单
        column(width = 4,
               card(
                 card_header("受试者清单 (支持多选对比)"),
                 div(
                   style = "height: 650px; overflow-y: auto;",
                   DTOutput("detail_table")
                 )
               )
        ),
        # 右侧 8格：Plotly 原生分面多线交互图
        column(width = 8,
               card(
                 card_header("选中受试者之全指标随访演变曲线 (多样本多色对比)"),
                 plotlyOutput("interactive_curve", height = "650px")
               )
        )
      )
    )
  )
)

# 4. Server 逻辑
server <- function(input, output, session) {

  # =======================================================
  # 面板 1：Sparklines 渲染
  # =======================================================
  output$summary_table <- renderReactable({
    render_advanced_spark <- function(values) {
      sparkline(
        unlist(values), type = "line", width = 110, height = 30,
        lineColor = "#2c3e50", fillColor = "#f4f6f7", lineWidth = 1.8,
        spotColor = "#3498db", minSpotColor = "#e74c3c", maxSpotColor = "#2ecc71", spotRadius = 2.5
      )
    }

    reactable(
      df_all_summary, pagination = TRUE, defaultPageSize = 5, pageSizeOptions = c(5, 10, 20), showPageSizeOptions = TRUE,
      compact = TRUE, bordered = TRUE, highlight = TRUE,
      defaultColDef = colDef(align = "center", headerStyle = list(background = "#f7f7f8", verticalAlign = "middle")),
      columns = list(
        subject_ID = colDef(name = "Subject ID", maxWidth = 120, style = list(fontWeight = "bold")),
        age = colDef(name = "Age", maxWidth = 100),
        group = colDef(name = "Group", maxWidth = 100)
      )
    )
  })

  # =======================================================
  # 面板 2：多维度联合筛选与 DT+Plotly 联动
  # =======================================================

  # 完美的联合响应式数据流：同时受 Age, Group 和 Red Mark 的控制
  filtered_raw_data <- reactive({
    req(input$selected_ages, input$selected_groups, input$selected_marks)
    df %>%
      filter(age %in% input$selected_ages) %>%
      filter(group %in% input$selected_groups) %>%
      filter(red_mark %in% input$selected_marks)
  })

  # 供左侧 DT 展示的去重唯一个体清单（保留 red_mark 用于表格展示）
  unique_subject_list <- reactive({
    req(filtered_raw_data())
    filtered_raw_data() %>%
      select(subject_ID, age, group, red_mark) %>%
      distinct(subject_ID, .keep_all = TRUE)
  })

  # 获取 DT 的当前选中的多个 Subject ID 向量
  selected_subject_ids <- reactive({
    sub_list <- unique_subject_list()
    if (is.null(sub_list) || nrow(sub_list) == 0) return(NULL)

    selected_idx <- input$detail_table_rows_selected

    # 门禁保护：如果没选，默认多选首位，保持右侧图形不崩溃
    if (is.null(selected_idx) || length(selected_idx) == 0) {
      return(sub_list$subject_ID[1])
    }

    valid_idx <- selected_idx[selected_idx <= nrow(sub_list)]
    if (length(valid_idx) == 0) return(sub_list$subject_ID[1])

    return(sub_list$subject_ID[valid_idx])
  })

  # 左侧唯一受试者列表：包含 Red Mark 展示
  output$detail_table <- renderDT({
    sub_list <- unique_subject_list()
    if (is.null(sub_list) || nrow(sub_list) == 0) return(NULL)

    # 重命名列用于前端 DT 优雅呈现
    display_df <- sub_list %>%
      rename(`Subject ID` = subject_ID, `Age` = age, `Group` = group, `Red Mark` = red_mark)

    datatable(
      display_df,
      selection = list(mode = "multiple", selected = 1), # 支持多选，默认选中第1行
      rownames = FALSE,
      options = list(
        dom = 't',
        deferRender = TRUE,
        scrollY = 580,
        paging = FALSE,
        ordering = FALSE,
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      )
    ) %>%
      # 选中行整行标蓝
      formatStyle(
        columns = 1:ncol(display_df),
        target = 'row',
        backgroundColor = styleEqual(c(TRUE, FALSE), c('#3498db', '#ffffff'), default = '#ffffff')
      )
  })

  # 右侧交互曲线：原生 plot_ly 驱动的多受试者连线对比图
  output$interactive_curve <- renderPlotly({
    target_ids <- selected_subject_ids()

    if (is.null(target_ids) || length(target_ids) == 0) {
      return(plotly_empty() %>% layout(title = list(text = "请在左侧选择受试者", x = 0.5)))
    }

    # 从全量原始数据中抓取所有被选中个体的全部随访记录进行绘图
    plot_data <- df %>%
      filter(subject_ID %in% target_ids) %>%
      select(subject_ID, visit, all_of(metrics_cols))

    if (nrow(plot_data) == 0) {
      return(plotly_empty() %>% layout(title = list(text = "所选受试者无合法记录", x = 0.5)))
    }

    # 转换为长面板数据
    plot_data <- plot_data %>%
      pivot_longer(cols = c(-subject_ID, -visit), names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = str_replace_all(Metric, "_", " ")) %>%
      filter(!is.na(Value))

    if (nrow(plot_data) == 0) return(NULL)

    # 用原生的 plot_ly 进行极速的多子图分面渲染
    subplot_list <- lapply(unique(plot_data$Metric), function(m) {
      sub_df <- plot_data %>% filter(Metric == m)

      plot_ly(sub_df, x = ~visit, y = ~Value, color = ~subject_ID, type = 'scatter', mode = 'lines+markers',
              text = ~paste0("ID: ", subject_ID, "<br>Visit: V", visit, "<br>Val: ", round(Value, 2)),
              hoverinfo = 'text', showlegend = (m == unique(plot_data$Metric)[1])) %>%
        layout(
          yaxis = list(title = m, titlefont = list(size = 10), tickfont = list(size = 8)),
          xaxis = list(title = "Visit", tickvals = 1:5, ticklabels = paste0("V", 1:5))
        )
    })

    # 合并为 3x3 交互式科研子图矩阵
    tryCatch({
      subplot(subplot_list, nrows = 3, shareX = TRUE, titleY = TRUE, margin = 0.05) %>%
        layout(
          title = list(text = "多指标个体演变并行对比面板", font = list(size = 14, weight = "bold")),
          legend = list(orientation = "h", x = 0, y = -0.1),
          margin = list(t = 60, b = 60, l = 40, r = 20)
        )
    }, error = function(e) {
      plotly_empty() %>% layout(title = list(text = "子图合并出错，请减少选中的受试者数量", x = 0.5))
    })

  })
}

shinyApp(ui, server)

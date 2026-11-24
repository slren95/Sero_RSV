library(shiny)
library(bslib)

# ============================================================
# 文件夹配置
# ============================================================

folder_paths <- c(
  "fit_main",
  "fit_sens1_0.05_0.05",
  "fit_sens1_0.10_0.10",
  "fit_sens1_0.10_0.15",
  "fit_sens1_2",
  "fit_sens2",
  "fit_sens3_3",
  "fit_sens3_9",
  "fit_sens3_12",
  "fit_sens4_5_5",
  "fit_sens4_1_1",
  "fit_sens5",
  "fit_main_new"
)


# 获取所有文件夹中的图片并集
get_all_images <- function(folders) {
  existing <- folders[dir.exists(folders)]
  if (length(existing) == 0) return(character(0))

  files <- unlist(lapply(existing, function(folder) {
    list.files(
      folder,
      pattern = "\\.(jpg|jpeg|png)$",
      ignore.case = TRUE,
      full.names = FALSE
    )
  }))

  sort(unique(files))
}


# ============================================================
# UI
# ============================================================

ui <- page_sidebar(
  title = "图片对比工具 - 多模型结果查看",
  theme = bs_theme(
    bootswatch = "flatly",
    primary = "#2c3e50"
  ),

  sidebar = sidebar(width = "30%",
    title = "控制面板",

    selectInput(
      "image_name",
      "选择图片",
      choices = NULL,
      width = "100%"
    ),

    hr(),

    checkboxGroupInput(
      "selected_folders",
      "选择要显示的模型",
      choices = setNames(folder_paths, folder_paths),
      selected = folder_paths[dir.exists(folder_paths)],
      width = "100%"
    ),

    hr(),

    div(
      style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
      h6("统计信息:", class = "text-muted"),
      textOutput("total_images"),
      textOutput("selected_count"),
      textOutput("current_image_status")
    ),

    hr(),

    actionButton(
      "refresh",
      "刷新文件列表",
      class = "btn-primary",
      style = "width: 100%;"
    )
  ),

  # 主界面：展示图片响应式网格
  uiOutput("image_grid")
)


# ============================================================
# SERVER
# ============================================================

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # ----------------------------------------------------------
  # 动态注册静态文件路径
  # ----------------------------------------------------------
  observe({
    for (folder in folder_paths) {
      if (dir.exists(folder)) {
        addResourcePath(
          prefix = basename(folder),
          directoryPath = normalizePath(folder)
        )
      }
    }
  })

  available_images <- reactiveVal(character(0))

  # ----------------------------------------------------------
  # 刷新逻辑（含调试打印）
  # ----------------------------------------------------------
  refresh_images <- function() {
    # 检查哪些文件夹真实存在
    existing <- folder_paths[dir.exists(folder_paths)]

    # 【调试输出】在 R 控制台打印检查结果
    message("--- 正在扫描文件夹 ---")
    message("当前工作目录: ", getwd())
    message("找到存在的文件夹数: ", length(existing))

    if (length(existing) == 0) {
      warning("错误：没有找到任何指定的文件夹！请检查工作目录或文件夹名称。")
      available_images(character(0))
      updateSelectInput(session, "image_name", choices = character(0))
      return()
    }

    # 读取所有图片（匹配 .jpg, .jpeg, .png, .bmp，不区分大小写）
    files <- unlist(lapply(existing, function(folder) {
      list.files(
        folder,
        pattern = "\\.(jpg|jpeg|png|bmp)$",
        ignore.case = TRUE, # 忽略大小写
        full.names = FALSE
      )
    }))

    images <- sort(unique(files))
    message("共扫描到图片数量: ", length(images))

    available_images(images)

    # 更新下拉框
    updateSelectInput(
      session,
      "image_name",
      choices = images,
      selected = if (length(images) > 0) images[1] else NULL
    )
  }

  # ----------------------------------------------------------
  # 保证 App 启动时【强制执行一次】刷新
  # ----------------------------------------------------------
  observeEvent(input$refresh, {
    refresh_images()
  }, ignoreInit = FALSE) # ignoreInit = FALSE 确保启动时就运行一次

  # 页面加载完成后主动刷新一次
  session$onFlushed(function() {
    refresh_images()
  }, once = TRUE)


  # ----------------------------------------------------------
  # 统计信息渲染
  # ----------------------------------------------------------
  output$total_images <- renderText({
    paste("总图片数:", length(available_images()))
  })

  output$selected_count <- renderText({
    paste("已选文件夹:", length(input$selected_folders), "/", length(folder_paths))
  })

  output$current_image_status <- renderText({
    if (is.null(input$image_name) || input$image_name == "") {
      return("未选择图片")
    }
    paste("当前查看:", input$image_name)
  })


  # ----------------------------------------------------------
  # 图片网格渲染
  # ----------------------------------------------------------
  output$image_grid <- renderUI({
    selected <- input$selected_folders

    req(length(selected) > 0)
    req(input$image_name)

    n <- length(selected)
    width_val <- switch(as.character(n),
                        "1" = 1,
                        "2" = 1/2,
                        "3" = 1/3,
                        1/2
    )

    card_list <- lapply(selected, function(folder) {
      img_name <- input$image_name
      img_path <- file.path(folder, img_name)

      if (file.exists(img_path)) {
        img_src <- paste0(basename(folder), "/", URLencode(img_name, reserved = TRUE))

        img_ui <- tags$img(
          src = img_src,
          style = "max-width: 100%; max-height: 550px; width: auto; height: auto; object-fit: contain;"
        )

        status <- tags$span(
          style = "font-size: 12px; color: #28a745; font-weight: bold;",
          "● 已找到"
        )
      } else {
        img_ui <- div(
          style = "width: 100%; height: 300px; display: flex; flex-direction: column; justify-content: center; align-items: center; background: #f8f9fa; border: 2px dashed #dee2e6; color: #6c757d;",
          div(style = "font-size: 36px; line-height: 1;", "∅"),
          div(style = "margin-top: 8px; font-size: 14px;", "该文件夹无此图片"),
          div(style = "font-size: 11px; margin-top: 4px; opacity: 0.7;", img_name)
        )

        status <- tags$span(
          style = "font-size: 12px; color: #dc3545; font-weight: bold;",
          "● 不存在"
        )
      }

      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          tags$strong(folder),
          status
        ),
        card_body(
          class = "d-flex align-items-center justify-content-center p-2",
          style = "min-height: 350px;",
          img_ui
        )
      )
    })

    layout_column_wrap(
      width = width_val,
      !!!card_list
    )
  })
}

# ============================================================
# RUN
# ============================================================

shinyApp(ui = ui, server = server)

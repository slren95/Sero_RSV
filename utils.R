
dir(path = 'code_02', pattern = '^code_02_.*[.]R$', full.names = TRUE) %>%
  walk(~{
    file.rename(.x,str_remove(.x,'truncat_extend_'))
  })

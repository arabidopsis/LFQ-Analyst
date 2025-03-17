progress_indicator <- function(text_message, session) {
  shiny::withProgress(
    message = text_message,
    detail = "Please wait for a while",
    value = 0,
    session = session,
    {
      for (i in 1:15) {
        shiny::incProgress(1 / 15)
        Sys.sleep(0.25)
      }
    }
  )
}

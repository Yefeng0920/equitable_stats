#-------------------------------------------------------------------------------------------------
# written by authors
#-------------------------------------------------------------------------------------------------
calc_evi <- function(m1, m2, sd1, sd2, n1, n2, rscale = 1) {
  # pooled standard deviation
  spooled <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  # Cohen's d
  d <- (m1 - m2) / spooled
  
  # Hedge's g
  J <- 1 - (3 / (4 * (n1 + n2 - 2) - 1))  # Correction factor
  g <- J * d
  
  # p-value
  sepooled <- sqrt((spooled^2 / n1) + (spooled^2 / n2))
  t_value <- (m1 - m2) / sepooled
  p_value <- 2 * pt(abs(t_value), df = n1 + n2 - 2, lower.tail = FALSE)
  
  # Bayes factor
  bayes_factor <- as.data.frame(meta.ttestBF(t = t_value, n1 = n1, n2 = n2, rscale = rscale))$bf
  
  # post-hoc power
  power <- pwr.t2n.test(n1 = n1, n2 = n2, d = d)$power
  
  # return a data frame with the results
  return(data.frame(
    d = d,
    g = g,
    p = p_value,
    bf = bayes_factor,
    pwr = power
  ))
}

#-------------------------------------------------------------------------------------------------
# modified from ihttps://github.com/davidsjoberg/ggsankey/blob/main/R/sankey.R
#-------------------------------------------------------------------------------------------------
utils::globalVariables(c(".", ".data", "x", "node", "next_node", "next_x", "..r"))
# importFrom(ggplot2, "%+replace%")
#' @importFrom ggplot2 %+replace%

# ** Support functions ----------
prepare_params <- function(...) {
  # Prepare aesthics for flow lines
  flow.aes <- list(...)
  removes <- names(flow.aes) %>%
    stringr::str_extract_all(., "(?<=flow.).*") %>% unlist()
  removes2 <- names(flow.aes) %>%
    stringr::str_subset(., "node") %>% unlist()
  flow.aes[c(removes, removes2)] <- NULL
  names(flow.aes) <- names(flow.aes) %>%
    stringr::str_replace_all("flow.", "")
  
  # Prepare aesthics for node boxes
  node.aes <- list(...)
  removes <- names(node.aes) %>%
    stringr::str_extract_all(., "(?<=node.).*") %>% unlist()
  removes2 <- names(node.aes) %>%
    stringr::str_subset(., "flow") %>% unlist()
  node.aes[c(removes, removes2)] <- NULL
  names(node.aes) <- names(node.aes) %>%
    stringr::str_replace_all(., "node.", "")
  
  return(list(flow.aes, node.aes))
}

find_default_space <- function(.df) {
  .df %>%
    dplyr::group_by(.data$n_x) %>%
    dplyr::summarise(n_groups = dplyr::n_distinct(.data$node),
                     freq = sum(.data$freq, na.rm = TRUE)) %>%
    dplyr::mutate(v = .data$freq / .data$n_groups / 4) %>%
    dplyr::pull(.data$v) %>%
    max()
}

sigmoid <- function(x_from, x_to, y_from, y_to, smooth = 5, n = 300) {
  x <- seq(-smooth, smooth, length = n)
  y <- exp(x) / (exp(x) + 1)
  out <- data.frame(x = (x + smooth) / (smooth * 2) * (x_to - x_from) + x_from,
                    y = y * (y_to - y_from) + y_from)
}


#long format -----------------------------------------------------------------
dlong <- function(.df, ..., value = NULL) {
  if("..r" %in% names(.df)) stop("The column name '..r' is not allowed")
  .vars <- dplyr::quos(...)
  
  if(!missing(value)) {
    value_var <- dplyr::enquo(value)
    out <- .df %>%
      dplyr::select(!!!.vars, value = !!value_var) %>%
      dplyr::mutate(..r = dplyr::row_number()) %>%
      tidyr::gather(x, node, -..r, -value) %>%
      dplyr::arrange(.data$..r) %>%
      dplyr::group_by(.data$..r) %>%
      dplyr::mutate(next_x = dplyr::lead(.data$x),
                    next_node = dplyr::lead(.data$node)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-..r) %>%
      dplyr::relocate(value, .after = dplyr::last_col())
  } else {
    out <- .df %>%
      dplyr::select(!!!.vars) %>%
      dplyr::mutate(..r = dplyr::row_number()) %>%
      tidyr::gather(x, node, -..r) %>%
      dplyr::arrange(.data$..r) %>%
      dplyr::group_by(.data$..r) %>%
      dplyr::mutate(next_x = dplyr::lead(.data$x),
                    next_node = dplyr::lead(.data$node)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-..r)
  }
  
  levels <- unique(out$x)
  
  out %>%
    dplyr::mutate(dplyr::across(c(x, next_x), ~factor(., levels = levels)))
}


#' @title sankey_themes
#' @name theme_sankey
#' @aliases theme_alluvial
#' @aliases theme_sankey_bump
#'
#' @description Minimal themes for sankey, alluvial and sankey bump plots
#'
#' @param base_size base font size, given in pts.
#' @param base_family base font family
#' @param base_line_size base size for line elements
#' @param base_rect_size base size for rect elements
#'
#' @export
theme_sankey <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black",
                                            size = ggplot2::rel(1)),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
  }

#' @rdname theme_sankey
#' @export
theme_alluvial <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
  }

#' @rdname theme_sankey
#' @export
theme_sankey_bump <-
  function(base_size = 11,
           base_family = "",
           base_line_size = base_size / 22,
           base_rect_size = base_size / 22) {
    {
      ggplot2::theme_bw(
        base_size = base_size,
        base_family = base_family,
        base_line_size = base_line_size,
        base_rect_size = base_rect_size
      ) %+replace%
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(
            fill = "white",
            colour = "transparent",
            size = ggplot2::rel(2)
          ),
          complete = TRUE,
          axis.line.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line("gray90")
        )
    }
  }


# FLOW LAYER ---------
StatSankeyFlow <- ggplot2::ggproto("StatSankeyFlow", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space", "smooth"),
                                   
                                   setup_data = function(data, params) {
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        flow_data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                          dplyr::summarise(flow_freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                        
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        flow_data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                          dplyr::summarise(flow_freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                        
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE),, .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      df <- data %>%
                                                        dplyr::left_join(flow_data, by = c("n_x", "node"))
                                                      
                                                      
                                                      
                                                      flows <- df %>%
                                                        dplyr::left_join(df %>%
                                                                           dplyr::select(n_x, node, ymin_end = ymin, ymax_end = ymax, xmin_end = xmin, xmax_end = xmax) %>%
                                                                           dplyr::distinct(),
                                                                         by = c("n_next_x" = "n_x", "next_node" = "node")) %>%
                                                        tidyr::drop_na(n_x, node, next_node, n_next_x, ymax_end, ymin_end, xmax_end, xmin_end) %>%
                                                        dplyr::mutate(r = dplyr::row_number()) %>%
                                                        dplyr::arrange(n_x, -r) %>%
                                                        dplyr::select(-r) %>%
                                                        dplyr::group_by(n_x, node) %>%
                                                        dplyr::mutate(cum_flow_freq = cumsum(flow_freq) - flow_freq) %>%
                                                        dplyr::ungroup() %>%
                                                        dplyr::group_by(n_x, n_next_x, node, next_node) %>%
                                                        dplyr::mutate(flow_start_ymax = ymax - cum_flow_freq,
                                                                      flow_start_ymin = flow_start_ymax - flow_freq)
                                                      
                                                      flows <- flows %>%
                                                        dplyr::arrange(n_x, n_next_x, next_node) %>%
                                                        dplyr::group_by(n_next_x, next_node) %>%
                                                        dplyr::mutate(cum_flow_freq_end = cumsum(flow_freq) - flow_freq) %>%
                                                        dplyr::mutate(flow_end_ymax = ymax_end - cum_flow_freq_end,
                                                                      flow_end_ymin = flow_end_ymax - flow_freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      flows <- flows %>%
                                                        dplyr::select(-n_x, -node, -freq, -ymax, -ymin, -xmin, -n_next_x, -next_node, -flow_freq, -ymin_end, -ymax_end, -xmax_end, -cum_flow_freq, -cum_flow_freq_end) %>%
                                                        dplyr::mutate(group = dplyr::row_number())
                                                      
                                                      flows %>%
                                                        dplyr::mutate(smooth = params$smooth) %>%
                                                        as.data.frame()
                                                    })
                                     
                                     
                                     
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     
                                     out1 <- sigmoid(data$xmax, data$xmin_end, data$flow_start_ymax, data$flow_end_ymax,
                                                     smooth = data$smooth)
                                     out2 <- sigmoid(data$xmin_end, data$xmax, data$flow_end_ymin, data$flow_start_ymin,
                                                     smooth = data$smooth)
                                     dplyr::bind_rows(out1, out2)
                                   }
)


# FLOW SANKEYBUMP LAYER ---------
StatSankeyBumpFlow <- ggplot2::ggproto("StatSankeyBumpFlow", ggplot2::Stat,
                                       extra_params = c("na.rm", "type", "space", "smooth"),
                                       
                                       setup_data = function(data, params) {
                                         
                                         purrr::map_dfr(unique(data$PANEL),
                                                        ~{
                                                          data <- data %>% dplyr::filter(PANEL == .x)
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(nodes = paste(node, x)) %>%
                                                            dplyr::arrange(x, -value) %>%
                                                            dplyr::mutate(bbb = dplyr::row_number()) %>%
                                                            dplyr::arrange(bbb) %>%
                                                            dplyr::mutate(nodes = fct_reorder(nodes, value, mean)) %>%
                                                            dplyr::arrange(node, x) %>%
                                                            dplyr::group_by(node) %>%
                                                            dplyr::mutate(next_x = dplyr::lead(x),
                                                                          node = nodes,
                                                                          next_node = dplyr::lead(nodes)) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::arrange(x, node)
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                          
                                                          if(!("value" %in% names(data))) {
                                                            flow_data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                              dplyr::summarise(flow_freq = dplyr::n(), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                            
                                                            data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::select(-n_next_x, -next_node) %>%
                                                              dplyr::group_by_all() %>%
                                                              dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                          } else {
                                                            flow_data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::group_by(n_x, node, n_next_x, next_node) %>%
                                                              dplyr::summarise(flow_freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                            
                                                            data <- data %>%
                                                              dplyr::mutate(group = 1) %>%
                                                              dplyr::select(-n_next_x, -next_node) %>%
                                                              dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                              dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                              dplyr::ungroup()
                                                          }
                                                          
                                                          if(is.null(params$space)) {
                                                            params$space <- find_default_space(data)
                                                          }
                                                          
                                                          data <- data %>%
                                                            dplyr::group_by(n_x) %>%
                                                            dplyr::arrange(node) %>%
                                                            dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                          ymin = ymax - freq) %>%
                                                            dplyr::ungroup()
                                                          
                                                          if(params$type == "sankey") {
                                                            data <- data %>%
                                                              dplyr::group_by(n_x) %>%
                                                              dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                            ymax = ymax - max(ymax)/2) %>%
                                                              dplyr::ungroup()
                                                          } else if (params$type == "alluvial"){
                                                            data <- data
                                                          }
                                                          
                                                          data <- data %>%
                                                            dplyr::mutate(xmin = n_x,
                                                                          xmax = n_x)
                                                          
                                                          df <- data %>%
                                                            dplyr::left_join(flow_data, by = c("n_x", "node"))
                                                          
                                                          flows <- df %>%
                                                            dplyr::left_join(df %>%
                                                                               dplyr::select(n_x, node, ymin_end = ymin, ymax_end = ymax, xmin_end = xmin, xmax_end = xmax, flow_freq_end = flow_freq) %>%
                                                                               dplyr::distinct(),
                                                                             by = c("n_next_x" = "n_x", "next_node" = "node")) %>%
                                                            tidyr::drop_na(n_x, node, next_node, n_next_x, ymax_end, ymin_end, xmax_end, xmin_end) %>%
                                                            dplyr::mutate(r = dplyr::row_number()) %>%
                                                            dplyr::arrange(n_x, -r) %>%
                                                            dplyr::select(-r) %>%
                                                            dplyr::group_by(n_x, node) %>%
                                                            dplyr::mutate(cum_flow_freq = cumsum(flow_freq) - flow_freq) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::group_by(n_x, n_next_x, node, next_node) %>%
                                                            dplyr::mutate(flow_start_ymax = ymax - cum_flow_freq,
                                                                          flow_start_ymin = flow_start_ymax - flow_freq)
                                                          
                                                          flows <- flows %>%
                                                            dplyr::arrange(n_x, n_next_x, next_node) %>%
                                                            dplyr::group_by(n_next_x, next_node) %>%
                                                            dplyr::mutate(cum_flow_freq_end = cumsum(flow_freq_end) - flow_freq_end) %>%
                                                            dplyr::mutate(flow_end_ymax = ymax_end - cum_flow_freq_end,
                                                                          flow_end_ymin = flow_end_ymax - flow_freq_end) %>%
                                                            dplyr::ungroup()
                                                          
                                                          flows <- flows %>%
                                                            dplyr::select(-n_x, -node, -freq, -ymax, -ymin, -xmin, -n_next_x, -next_node, -flow_freq, -ymin_end, -ymax_end, -xmax_end, -cum_flow_freq, -cum_flow_freq_end) %>%
                                                            dplyr::mutate(group = dplyr::row_number())
                                                          
                                                          flows %>%
                                                            rowwise() %>%
                                                            dplyr::mutate(..groupqq = stringr::str_remove(nodes, as.character(x))) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::group_by(..groupqq) %>%
                                                            dplyr::mutate(group = dplyr::cur_group_id()) %>%
                                                            dplyr::ungroup() %>%
                                                            dplyr::select(-..groupqq) %>%
                                                            dplyr::mutate(smooth = params$smooth) %>%
                                                            as.data.frame()
                                                        })
                                       },
                                       
                                       compute_group = function(data, scales) {
                                         
                                         out1 <- purrr::map_dfr(1:nrow(data), ~{
                                           datat <- data %>% dplyr::slice(.x)
                                           sigmoid(datat$xmax, datat$xmin_end, datat$flow_start_ymax, datat$flow_end_ymax,
                                                   smooth = datat$smooth)
                                         }) %>%
                                           dplyr::arrange(x)
                                         out2 <- purrr::map_dfr(1:nrow(data), ~{
                                           datat <- data %>% dplyr::slice(.x)
                                           sigmoid(datat$xmin_end, datat$xmax, datat$flow_end_ymin, datat$flow_start_ymin,
                                                   smooth = datat$smooth)
                                         }) %>%
                                           dplyr::arrange(-x)
                                         
                                         dplyr::bind_rows(out1, out2)
                                       }
)

# TEXT LAYER -------
StatSankeyText <- ggplot2::ggproto("StatSankeyText", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space"),
                                   
                                   setup_data = function(data, params) {
                                     
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(x = n_x,
                                                                      y = ymin + (ymax - ymin)/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      
                                                      return(as.data.frame(data))
                                                    })
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     data
                                   }
)


# NODE LAYER -------
StatSankeyNode <- ggplot2::ggproto("StatSankeyNode", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space", "smooth"),
                                   
                                   setup_data = function(data, params) {
                                     
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      
                                                      data <- data %>% dplyr::filter(PANEL == .x)
                                                      data <- data %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))
                                                      
                                                      if(!("value" %in% names(data))) {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_all() %>%
                                                          dplyr::summarise(freq = dplyr::n(), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      } else {
                                                        data <- data %>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                          dplyr::group_by_at(dplyr::vars(dplyr::everything(), -value)) %>%
                                                          dplyr::summarise(freq = sum(value, na.rm = TRUE), .groups = "keep") %>%
                                                          dplyr::ungroup()
                                                      }
                                                      
                                                      if(is.null(params$space)) {
                                                        params$space <- find_default_space(data)
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::group_by(n_x) %>%
                                                        dplyr::mutate(ymax = cumsum(freq) + (dplyr::row_number() - 1)*params$space,
                                                                      ymin = ymax - freq) %>%
                                                        dplyr::ungroup()
                                                      
                                                      if(params$type == "sankey") {
                                                        data <- data %>%
                                                          dplyr::group_by(n_x) %>%
                                                          dplyr::mutate(ymin = ymin - max(ymax)/2,
                                                                        ymax = ymax - max(ymax)/2) %>%
                                                          dplyr::ungroup()
                                                      } else if (params$type == "alluvial"){
                                                        data <- data
                                                      }
                                                      
                                                      data <- data %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      
                                                      if("shift" %in% names(data)) {
                                                        data <- data %>%
                                                          dplyr::mutate(dplyr::across(dplyr::contains("y"), ~ . + shift))
                                                      }
                                                      
                                                      return(as.data.frame(data))
                                                    })
                                     
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     data
                                   }
)


sankey_p <- function(mapping = NULL,
                     data = NULL,
                     position = "identity",
                     na.rm = FALSE,
                     show.legend = NA,
                     space = NULL,
                     type = "sankey",
                     width = .1,
                     smooth = 8,
                     inherit.aes = TRUE,
                     ...
) {
  params_list <- prepare_params(...)
  
  list(
    flow = ggplot2::layer(
      stat = StatSankeyFlow,
      data = data,
      mapping = mapping,
      geom = "polygon",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          smooth = smooth,
          type = type,
          params_list[[1]]
        )
      )
    ),
    
    node = ggplot2::layer(
      stat = StatSankeyNode,
      data = data,
      mapping = mapping,
      geom = ggplot2::GeomRect,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          smooth = smooth,
          type = type,
          params_list[[2]]
        )
      )
    )
  )
  
  
}


sankey_p_label <- function(mapping = NULL,
                           data = NULL,
                           position = "identity",
                           na.rm = FALSE,
                           show.legend = NA,
                           space = NULL,
                           type = "sankey",
                           width = .1,
                           inherit.aes = TRUE,
                           ...) {
  # Prepare aesthics for label
  label.aes <- list(...)
  
  list(
    label = ggplot2::layer(
      stat = StatSankeyText,
      data = data,
      mapping = mapping,
      geom = "label",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          type = type,
          label.aes
        )
      )
    )
  )
}
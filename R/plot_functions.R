#' Trend plot
#' @param y_0 True predict target, length n vector
#' @param y_hat Predicted values by m different methods, n-by-m matrix, colnames well defined
#' @param dates Date vector of class "Date". \cr
#' Use \textit{lubridate} and \textit{zoo} to transfer original date vector to "Date" class
#' @param pt_num Number date points (pt_num + 1) wanted on the x-axis
#' @param col_vec A vector indicates colors of different trends. \cr
#' e.g.:  c("#D8DBE2", "#F46B7B", "#518DE8", "#FFBC42" )\cr
#' If NULL, use ggplot default.
#' @param line_size line size for each trend
#' @param alpha_size degree of appearance
#' @param xlab x-axis label
#' @param ylab y-axis label
#'
#' @return A ggplot output
#'
#' @examples
#' data("forecast_result")
#' data("raw_data_h1")
#' D <- raw_data_h1[!is.na(raw_data_h1$LongReturn), ]
#' dates <- zoo::as.Date.yearmon(D$yyyymm[-(1:180)])
#' y0 <- forecast_result$y[-(1:180)]
#' y_hat <-cbind(forecast_result$Lasso$y_hat,
#'               forecast_result$Lasso_Std$y_hat,
#'               forecast_result$ALasso$y_hat)
#' colnames(y_hat) <- c("lasso", "lasso_std", "alasso")
#' plot_trend(y0, y_hat, dates)
#'
#' @export plot_trend


plot_trend <- function(y_0,
                       y_hat,
                       dates,
                       pt_num = 4,
                       col_vec = NULL,
                       line_size = rep(0.75, ncol(y_hat) + 1),
                       alpha_size = rep(0.75, ncol(y_hat) + 1),
                       xlab = NULL,
                       ylab = NULL) {


    n <- length(y_0)
    m <- ncol(y_hat)
    methods_use <- colnames(y_hat)

    if(class(dates) != "Date"){
        dates <- zoo::as.Date(dates)
    }

    if(is.null(col_vec)){
        col_vec <- gg_color_hue(m + 1)
    }

    # Prepare data frame
    df_plot <- cbind(data.frame(date = dates, y_0 = y_0), as.data.frame(y_hat))
    df_plot <- reshape2::melt(df_plot, id = "date")

    pt <- dates[seq(from = 1, to = n, by = (n / pt_num - 1))]

    p_out <- ggplot(data = df_plot) +
        geom_line(mapping = aes(
            x = date,
            y = value,
            color = variable,
            size = variable,
            alpha = variable
        )) +
        scale_colour_manual(breaks = c("y_0", methods_use),
                            values = col_vec,
                            labels = c("True value", methods_use),
                            guide = guide_legend(override.aes = aes(alpha = NA))) +
        scale_size_manual(breaks = c("y_0", methods_use),
                          values = line_size,
                          labels = c("True value", methods_use)) +
        scale_alpha_manual(breaks = c("y_0", methods_use),
                           values = alpha_size,
                           labels = c("True value", methods_use)) +
        scale_x_date(breaks = pt, labels = as.character(lubridate::year(pt))) +
        labs(x = xlab, y = ylab) +
        theme(
            panel.background =  element_blank(),
            panel.border = element_rect(
                linetype = 1,
                colour = "black",
                fill = NA
            ),
            panel.grid.major = element_line(linetype = 2, color = "grey90"),
            legend.title = element_blank(),
            legend.position = "bottom"
        )

    return(p_out)

}

#' Coefficient plot
#'
#' @param coef_est estimated slope by m different methods, list of length m with each element n-by-p matrix\cr
#' names of list well defined
#' @param dates Date vector of class "Date". \cr
#' @param pt_num Number date points (pt_num + 1) wanted on the x-axis
#' Use \textit{lubridate} and \textit{zoo} to transfer original date vector to "Date" class
#' @param col_vec A vector indicates colors of different trends. \cr
#' e.g.:  c("#D8DBE2", "#F46B7B", "#518DE8", "#FFBC42" )\cr
#' If NULL, use ggplot default.
#' @param line_size line size for each trend
#' @param alpha_size degree of appearance
#' @param xlab x-axis label
#' @param ylab y-axis label
#'
#' @return A ggplot output
#'
#' @examples
#' data("forecast_result")
#' data("raw_data_h1")
#' D <- raw_data_h1[!is.na(raw_data_h1$LongReturn), ]
#' dates <- zoo::as.Date.yearmon(D$yyyymm[-(1:180)])
#'
#' coef_est <- list(lasso = forecast_result$Lasso$beta_hat[, 1:6],
#'                  alasso = forecast_result$ALasso$beta_hat[, 1:6])
#' plot_coef(coef_est, dates)
#'
#' @export plot_coef
#'
plot_coef <- function(coef_est, dates,
                      pt_num = 4,
                      col_vec = NULL,
                      line_size = rep(0.75, length(coef_est)),
                      alpha_size = rep(0.75, length(coef_est)),
                      xlab = NULL,
                      ylab = NULL,
                      num_col = 1){

    n <- nrow(coef_est[[1]])
    p <- ncol(coef_est[[1]])
    m <- length(coef_est)

    methods_use <- names(coef_est)
    var_names <- colnames(coef_est[[1]])

    if(class(dates) != "Date"){
        dates <- zoo::as.Date(dates)
    }

    if(is.null(col_vec)){
        col_vec <- gg_color_hue(m)
    }

    pt <- dates[seq(from = 1, to = n, by = (n / pt_num - 1))]

    df_plot <- NULL
    for(i in 1:m){

        method_i <- methods_use[i]
        df_temp <- cbind(data.frame(date = dates), as.data.frame(coef_est[[method_i]]))
        df_temp$method <- method_i

        df_plot <- rbind(df_plot, df_temp)
    }

    df_plot <- reshape2::melt(df_plot, id = c("date","method"))
    df_plot$variable = factor(df_plot$variable, levels = var_names)

    p_out <- ggplot(data = df_plot) +
        geom_line(mapping = aes(x = date, y = value, color = method, size = method, alpha = method)) +
        scale_x_date(breaks = pt, labels = as.character(lubridate::year(pt))) +
        labs(x = xlab, y = ylab) +
        theme(
            strip.background = element_blank(),
            strip.text =  element_text(face = "bold"),
            panel.background =  element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_rect(
                linetype = 1,
                colour = "black",
                fill = NA
            ),
            panel.grid.major = element_line(linetype = 2, color = "grey90"),
            legend.position = "bottom",
            legend.title = element_blank()
        ) +
        scale_color_manual(values = col_vec,
                           guide = guide_legend(override.aes=aes(alpha=NA))) +
        scale_size_manual(values = line_size) +
        scale_alpha_manual(values = alpha_size) +
        facet_wrap(~variable, ncol = num_col, scales = "free", strip.position="right")

    return(p_out)

}

#' Share a legend between multiple plots using grid.arrange
#' from https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#'
#' @export grid_arrange_shared_legend

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {


    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

    grid.newpage()
    grid.draw(combined)

    # return gtable invisibly
    invisible(combined)

}




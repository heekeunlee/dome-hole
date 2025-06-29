# 자동 혼합 최적화 타공 시뮬레이터 (2~5개 조합, AI 추천 포함)

if (!require("shiny")) install.packages("shiny")
if (!require("plotly")) install.packages("plotly")

library(shiny)
library(plotly)

hole_shapes <- c("circle", "square", "star", "flower", "ellipse", "cross", "trefoil")

# 기본 곡면 생성 함수
generate_edge_curve <- function(x, y, L, W, T, H) {
  cx <- L / 2; cy <- W / 2
  dx <- abs(x - cx); dy <- abs(y - cy)
  rx <- (L - T) / 2; ry <- (W - T) / 2
  over_x <- pmax(dx - rx, 0); over_y <- pmax(dy - ry, 0)
  edge <- pmin(1, sqrt(over_x^2 + over_y^2) / (T / 2))
  z <- H * (1 - edge^4); z[z < 0] <- 0; return(z)
}

make_hole_mask <- function(x_seq, y_seq, hx, hy, r, shape_fn) {
  outer(x_seq, y_seq, function(x, y) shape_fn(x, y, hx, hy, r))
}

get_shape_function <- function(shape) {
  switch(shape,
         "circle" = function(x, y, hx, hy, r) (x - hx)^2 + (y - hy)^2 <= r^2,
         "square" = function(x, y, hx, hy, r) abs(x - hx) <= r & abs(y - hy) <= r,
         "star" = function(x, y, hx, hy, r) (abs(x - hx)^1.5 + abs(y - hy)^1.5) <= r^1.5,
         "flower" = function(x, y, hx, hy, r) {
           theta <- atan2(y - hy, x - hx)
           r_eff <- r * (0.6 + 0.4 * sin(5 * theta))
           sqrt((x - hx)^2 + (y - hy)^2) <= r_eff
         },
         "ellipse" = function(x, y, hx, hy, r) ((x - hx)^2 / r^2 + (y - hy)^2 / (r * 0.6)^2) <= 1,
         "cross" = function(x, y, hx, hy, r) (abs(x - hx) <= r * 0.2 | abs(y - hy) <= r * 0.2),
         "trefoil" = function(x, y, hx, hy, r) {
           theta <- atan2(y - hy, x - hx)
           r_eff <- r * (1 + 0.3 * sin(3 * theta))
           sqrt((x - hx)^2 + (y - hy)^2) <= r_eff
         })
}

compute_stress <- function(z) {
  zx <- apply(z, 2, function(col) c(NA, diff(col, differences = 2), NA))
  zy <- apply(z, 1, function(row) c(NA, diff(row, differences = 2), NA))
  zx <- t(zx)
  stress <- sqrt(zx^2 + zy^2)
  stress[is.na(stress)] <- 0
  stress / max(stress, na.rm = TRUE)
}

generate_combinations <- function(shapes) {
  n <- length(shapes)
  all_combinations <- expand.grid(rep(list(seq(0, 1, by = 0.1)), n))
  all_combinations <- all_combinations[rowSums(all_combinations) == 1, ]
  colnames(all_combinations) <- shapes
  return(all_combinations)
}

ui <- fluidPage(
  titlePanel("혼합 타공 최적화 시뮬레이터"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("L", "패널 길이 (mm)", 30, 60, 40),
      sliderInput("W", "패널 너비 (mm)", 30, 60, 40),
      sliderInput("T", "가장자리 곡률폭 (mm)", 5, 30, 10),
      sliderInput("H", "중앙 높이 (mm)", 1, 5, 2.5),
      sliderInput("D", "타공 지름 (mm)", 0.2, 1, 0.5, step = 0.1),
      sliderInput("S", "타공 간격 (mm)", 0.5, 3, 1, step = 0.1),
      checkboxGroupInput("selected_shapes", "혼합할 타공 모양 (2~5개 선택)",
                         choices = hole_shapes, selected = c("circle", "flower")),
      actionButton("run", "최적 혼합 계산 실행")
    ),
    mainPanel(
      plotlyOutput("heatmap"),
      textOutput("best_mix")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$run, {
    req(length(input$selected_shapes) >= 2)
    x_seq <- seq(0, input$L, length.out = 100)
    y_seq <- seq(0, input$W, length.out = 100)
    z_mat <- outer(x_seq, y_seq, Vectorize(function(x, y) {
      generate_edge_curve(x, y, input$L, input$W, input$T, input$H)
    }))
    base_stress <- compute_stress(z_mat)
    base_mean <- mean(base_stress, na.rm = TRUE)
    
    hole_x <- seq(input$D / 2, max(x_seq), by = input$S)
    hole_y <- seq(input$D / 2, max(y_seq), by = input$S)
    holes <- expand.grid(hx = hole_x, hy = hole_y)
    combinations <- generate_combinations(input$selected_shapes)
    
    reliefs <- numeric(nrow(combinations))
    stress_list <- list()
    for (i in 1:nrow(combinations)) {
      combo <- combinations[i, ]
      mask_total <- matrix(1, nrow = length(x_seq), ncol = length(y_seq))
      shape_ids <- sample(rep(names(combo), round(combo * nrow(holes))))
      shape_ids <- shape_ids[1:nrow(holes)]
      for (j in 1:nrow(holes)) {
        shape_fn <- get_shape_function(shape_ids[j])
        mask <- make_hole_mask(x_seq, y_seq, holes$hx[j], holes$hy[j], input$D / 2, shape_fn)
        mask_total[mask] <- mask_total[mask] * 0.5
      }
      stress_mod <- base_stress * mask_total
      reliefs[i] <- (base_mean - mean(stress_mod, na.rm = TRUE)) / base_mean
      stress_list[[i]] <- stress_mod
    }
    best_i <- which.max(reliefs)
    best_combo <- combinations[best_i, ]
    best_stress <- stress_list[[best_i]]
    
    x_zoom_range <- c(input$L - 15, input$L)
    y_zoom_range <- c(0, 10)
    over_75_ratio <- round(sum(best_stress >= 0.75, na.rm = TRUE) / length(best_stress) * 100, 2)
    avg_stress_ratio <- round(mean(best_stress, na.rm = TRUE) * 100, 2)
    
    output$heatmap <- renderPlotly({
      plot_ly(
        z = best_stress,
        x = x_seq,
        y = y_seq,
        type = "heatmap",
        colorscale = list(c(0, "blue"), c(0.25, "green"), c(0.5, "yellow"), c(0.75, "orange"), c(1, "red")),
        zmin = 0, zmax = 1,
        showscale = TRUE
      ) %>%
        add_contour(
          z = best_stress,
          x = x_seq,
          y = y_seq,
          contours = list(start = 0, end = 1, size = 0.05),
          line = list(width = 0.3, color = "gray"),
          showscale = FALSE
        ) %>%
        layout(
          title = "우측 하단 모서리 응력 분포 확대",
          xaxis = list(title = "X (mm)", range = x_zoom_range, scaleanchor = "y"),
          yaxis = list(title = "Y (mm)", range = y_zoom_range),
          annotations = list(
            list(
              x = 1.08, y = 0.5, xref = "paper", yref = "paper",
              text = paste0("응력 ≥ 75%: ", over_75_ratio, "%"),
              showarrow = FALSE, font = list(size = 12)
            ),
            list(
              x = 1.08, y = 0.4, xref = "paper", yref = "paper",
              text = paste0("평균 응력: ", avg_stress_ratio, "%"),
              showarrow = FALSE, font = list(size = 12)
            )
          )
        )
    })
    
    output$best_mix <- renderText({
      paste("최적 혼합 비율:", paste(names(best_combo), round(100 * best_combo), "%", collapse = ", "))
    })
  })
}

shinyApp(ui, server)
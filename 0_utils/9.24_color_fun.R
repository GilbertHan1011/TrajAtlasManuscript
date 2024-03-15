
#== tutorial at https://nanx.me/blog/post/ggplot2-color-interpolation/
pal_ramp <- function(values) {
  force(values)
  function(n) {
    if (n <= length(values)) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}
pal_adaptive <- function(name, palette, alpha = 1) {
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db[[name]][[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  pal_ramp(unname(alpha_cols))
}

scale_color_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("colour", name, pal_adaptive(name, palette, alpha), ...)
}

scale_fill_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("fill", name, pal_adaptive(name, palette, alpha), ...)
}

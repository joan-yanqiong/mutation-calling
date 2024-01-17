#' Create venndiagram
#' @inheritParams  VennDiagram::venn.diagram
#' @return NA
#' @export
#' @importFrom VennDiagram venn.diagram
#' @importFrom glue glue
#' @importFrom colorjam rainbowJam
create_venndiagram <- function(x, category.names, filename, main = "", main.cex = 0.3, cat.cex = 0.15, height = 480, width = 600, ...) {
    venn.diagram(
        x = x,
        category.names = category.names,
        filename = filename,
        output = FALSE,
        main = main,
        disable.logging = FALSE,
        # Output features
        imagetype = "png",
        height = height,
        width = width,
        resolution = 600,
        compression = "lzw",

        # Circles
        lwd = 2,
        lty = "blank",
        fill = rainbowJam(length(category.names)),

        # Numbers
        cex = .4,
        fontface = "bold",
        fontfamily = "sans",

        # Set names
        cat.cex = cat.cex,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.fontfamily = "sans",

        # Title
        main.fontfamily = "sans",
        main.cex = main.cex,
        main.fontface = "bold",
    )
}

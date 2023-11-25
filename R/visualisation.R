#' Create venndiagram
#' @inheritParams  VennDiagram::venn.diagram
#' @return NA
#' @export
#' @importFrom VennDiagram venn.diagram
#' @importFrom glue glue
#' @importFrom colorjam rainbowJam
create_venndiagram <- function(x, category.names, filename, main = "", main.cex = 0.3, cat.cex = 0.15, ...) {
    venn.diagram(
        x = x,
        category.names = category.names,
        filename = filename,
        output = FALSE,
        main = main,
        disable.logging = FALSE,
        # Output features
        imagetype = "png",
        height = 480,
        width = 480,
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
        main.fontface = "bold"
    )
}

# dagitty, tidygraph converters

dag_to_tidy_graph <- function(dag){
    # Extract edges
    edges <- as.data.frame( dagitty::edges(dag) )
    edges$edge_w = 1
    nodes <- tibble(name = unique(c(edges$v, edges$w))) 
    
    # Build tidygraph
    tg <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
    return(tg)
}

dagitty_from_edge_df <- function(edges) {
    # edges is a df with "from" and "to"
    dagitty(str_c(
        "dag{", 
        str_c(edges |> 
                  mutate(edge_string = str_c(from, "->", to)) |> 
                  pull(edge_string),
              collapse = ";\n")
        ,"}"))
}

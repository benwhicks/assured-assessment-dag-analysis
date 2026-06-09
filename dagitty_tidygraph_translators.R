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

tidy_graph_to_dag <- function(g) {
    nodes <- g %N>% as_tibble() 
    edges <- g %E>% as_tibble() 
    
    edge_strings <- 
        edges |> mutate(
        from.n = nodes$name[from],
        to.n = nodes$name[to]
        ) |> 
        mutate(edge_string = str_c(from.n, "->", to.n, sep = " ")) |> 
        pull(edge_string)
    
    dag_string <- str_c("dag{", 
                        str_c(edge_strings, collapse = "\n"),
                        "}")
    
    return(dagitty::dagitty(dag_string))
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

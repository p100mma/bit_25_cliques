critical_mst_thr<- function(S) {
#igraph gives minimum spanning tree so we switch to dissimilarity graph
dS= max(S) - S
G_ds= igraph::graph_from_adjacency_matrix(adjmatrix=dS, mode="undirected", weighted=TRUE, diag=FALSE)

igraph::mst(G_ds)-> ds_mst
#convert back to similarities so we sort weights of MAXIMIMUM spanning tree as in description
sort( max(S) - igraph::E(ds_mst)$weight, decreasing=TRUE)-> wMSTsorted
wMSTsorted<- wMSTsorted[ wMSTsorted!=0 ]
wdiffs= -1*(diff(wMSTsorted))
wmeans= (wMSTsorted[ 2: length(wMSTsorted) ] + wMSTsorted[ 1: (length(wMSTsorted)-1) ])/2
wmeans[[ which.max(wdiffs) ]]-> w_t
attr(w_t, "mst_weights")<- wMSTsorted
attr(w_t, "wdiffs")<- wdiffs
attr(w_t, "wmeans")<- wmeans
w_t
}

init_potentialPerClique<- function(S, P){
 uql=unique(P)
 uql<- uql[ uql!=0 ]
 labs=sort(uql,decreasing=FALSE) 
 stopifnot(all(labs==1:length(uql)))
 potential_indicator<- matrix( nrow=length(uql),
			       ncol=length(P))
 for (i in labs){
	potential_indicator[i,]=colAlls(S[P==i,!(P %in% c(0,i))]>0)
 }
 potential_indicator
}

initialize_state<- function(S, P){
	state<- environment()
	state$P<-P
	state$S-S
	state$potential_indicator<- init_potentialPerClique(S,P)
	return(state)
}

link_strength<- function(v_connections, clique_number, P, aggregation_fun=mean){
	aggregation_fun(v_connections[P==clique_number])
}

calculate_preferences<- function(prev_state){

preferences<-lapply(prev_state$potential_indicator, function(v_mask) if (all(!v_mask)) return(NA)


}


join_by_strongest<- function(prev_state){



}




# Help functions
if(FALSE){
create.neighbourhood.array <- function(m = c(5, 5), sp = 2, type = "queen", torus = TRUE){
  # Checking input:
  if(length(m) == 1) m <- c(m, 1)
  if(!is.numeric(m)) stop("m must be numeric.")
  if(!is.numeric(sp)) stop("sp must be numeric.")
  if(!(type %in% c("rook","queen"))) stop("Type must be either rook or queen")
  if(!is.logical(torus)) stop("Torus must be logical.")
  
  # Neighbourmatrix:
  W <- spdep::cell2nb(m[1],m[2], type = type, torus = torus)
  # Neighbour array:
  Warr <- array(0, dim = c(prod(m),prod(m), sp+1))
  Warr[,,1] <- diag(prod(m))
  if(sp == 0) return(Warr)
  W <- spdep::nblag(W, max(sp,2))
  for(i in 1:(sp)){
    Warr[,,i+1] <- spdep::nb2mat(W[[i]])
    Warr[,,i+1] <- Warr[,,i+1]/max(Warr[,,i+1])
  }
  #Warr[,,2]<- Warr[,,2]-diag(prod(m))
  return(Warr)
}
}
correct.names <- function(names){
  if(all(table(names)==1)) return(names)
  k <- table(names)
  k <- k[which(k>1)]
  for(i in names(k)){
    names[which(names==i)] <- paste(i, 1:k[i], sep="")
  }
  names
}


read_cell <- function(filename){
  df1<-as.matrix(read.table(filename, header=FALSE, sep=",", stringsAsFactors = FALSE))
  rownames(df1)<-1:nrow(df1)
  colnames(df1)<-1:ncol(df1)
  return(df1)
}

plotting_df <- function(df, main = NULL, path = NULL, hjust=0.5, midpoint = mean(df), limits = NULL, guides = TRUE){
  #df <- abs(df - mean(df[!is.na(df)]))
  #df <- t(apply(df, 1, diff))
  #df <- log(df)
  df.melt <- melt(df)
  #guide <- ifelse(guides, , FALSE)
  if(is.null(limits))limits <- c(-1,1)*max(abs(df.melt$value))
  (p <- ggplot(df.melt, aes(Var2*10, Var1))+geom_raster(aes(fill=value))+
      scale_fill_gradient2(low="blue", mid = "white", high ="red", limits = limits,
                           midpoint = midpoint,
                           name = "", guide = guide_colourbar(barheight = 15, ticks.colour = "grey40"))+
      theme_bw()+xlab("Time (sec)")+ylab("Space")+ggtitle(main)+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust=hjust),
            #axis.line = element_line(colour = "black"),
            panel.border = element_blank())+
      scale_x_discrete(expand =c(0,0), limits = seq(0,400*10+100,100*10))+
      scale_y_discrete(expand =c(0,0), limits = seq(0,nrow(df), 20)))
  if(!is.null(path))
    ggsave(path, width = 9, height =4)
  return(p)
}

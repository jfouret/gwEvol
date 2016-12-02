venn2=function(names=c("",""),col=c('red','blue'),fontsize=16,weigths=c("","","")){
  #16 recommended for A5 saving
  require(ggplot2)
  require(gridExtra)
  require(ggthemes)
  require(cowplot)
  cercle=function(col='black',alpha=0.5){
    y=1
    data=as.data.frame(y)
    data$class="venn"
    ggplot(data,aes(x='',y=y))+
      coord_polar()+
      geom_bar(width = 1,stat='identity',fill=col,alpha=alpha)+
      theme_void()
  }
  u=1/4
  venn=ggdraw() +
    draw_plot(cercle(col[1]), 0, 0, 3*u, 1) +
    draw_plot(cercle(col[2]), u, 0, 3*u, 1)+
    draw_text(text = names[1],x = 1.5*u,y = 1,size = fontsize,vjust=1.5)+
    draw_text(text = names[2],x = 2.5*u,y = 1,size = fontsize,vjust=1.5)+
    draw_text(text = names[3],x = 2.5*u,y = 1,size = fontsize,vjust=1.5)+
    draw_text(text = weigths[1],x = u,y = 0.5,size = 2*fontsize)+
    draw_text(text = weigths[2],x = 2*u,y = 0.5,size = 2*fontsize)+
    draw_text(text = weigths[3],x = 3*u,y = 0.5,size = 2*fontsize)
  return(venn)
}

venn3=function(names=c("","",""),col=c('red','blue','yellow'),fontsize=16,weigths=c("","","",'','','','')){
  #16 recommended for 5x5 inches saving
  require(ggplot2)
  require(gridExtra)
  require(ggthemes)
  require(cowplot)
  cercle=function(col='black',alpha=0.5){
    y=1
    data=as.data.frame(y)
    data$class="venn"
    ggplot(data,aes(x='',y=y))+
      coord_polar()+
      geom_bar(width = 1,stat='identity',fill=col,alpha=alpha)+
      theme_void()
  }
  u=1/4
  venn=ggdraw() +
    draw_plot(cercle(col[1]), 0, u/2, 3*u, 1)+
    draw_plot(cercle(col[2]), u, u/2, 3*u, 1)+
    draw_plot(cercle(col[3]), u/2, -u/2, 3*u, 1)+
    draw_text(text = names[1],x = 1.5*u,y = 1,size = fontsize,vjust=1.5)+
    draw_text(text = names[2],x = 2.5*u,y = 1,size = fontsize,vjust=1.5)+
    draw_text(text = names[3],x = 2*u,y = 0,size = fontsize,vjust=-1)+
    draw_text(text = weigths[1],x = u,y = u*3,size = 2*fontsize)+
    draw_text(text = weigths[3],x = 2*u,y = u*0.75,size = 2*fontsize)+
    draw_text(text = weigths[2],x = 3*u,y = u*3,size = 2*fontsize)+
    draw_text(text = weigths[6],x = 1.25*u,y = u*1.75,size = 2*fontsize)+
    draw_text(text = weigths[5],x = 2.75*u,y = u*1.75,size = 2*fontsize)+
    draw_text(text = weigths[4],x = 2*u,y = u*3.1,size = 2*fontsize)+
    draw_text(text = weigths[7],x = 2*u,y = u*2.1,size = 2*fontsize)
  return(venn)
}



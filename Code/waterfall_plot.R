# Code for waterfall plot

waterfall_f <- function(df, offset=0.3) {
  
  library(ggplot2)
  library(scales)
  library(dplyr)
  
  ## Add the order column to the raw data frame and order appropriately
  #df <- df %>% mutate(order=as.numeric(category)) %>% arrange(order)
  df$order <- 1:length(df$category)
  
  ## The last value needs to be negated so that it goes down to
  ## zero.  Throws a warning if the cumulative sum doesn't match.
  last.id <- nrow(df)
  df$value[last.id] <- -df$value[last.id]
  
  ## Calculate the cumulative sums
  df <- df %>% mutate(TZA_2010=cumsum(value))
  
  ## Throw a warning if the values don't match zero as expected
  final_value <- tail(df$TZA_2010, 1)
  #   if (final_value!=0) {
  #     warning(sprintf("Final value doesn't return to 0.  %.2d instead.", final_value))
  #   }
  
  ## Calculate the max and mins for each category and sector
  df <- transform(df, min.val=c(0, head(TZA_2010, -1)),
                  max.val=c(head(TZA_2010, -1), 0))    
  df <- df %>% group_by(order, category, sector, value, TZA_2010) %>%
    summarize(min=min(min.val, max.val), max=max(min.val, max.val))
  
  ## Create the lines data frame to link the bars
  lines <- df %>% group_by(order) %>% summarize(cs=max(TZA_2010))
  lines <- with(lines, data.frame(x=head(order, -1),
                                  xend=tail(order, -1),
                                  y=head(cs, -1),
                                  yend=head(cs, -1)))
  
  
  ## Add the offset parameter
  df <- transform(df, offset=offset)
  
  ## Make the plot    
  gg <- ggplot() +
    geom_segment(data=lines, aes(x=x, y=y, xend=xend, yend=yend), linetype="dashed")  +
    geom_rect(data=df, aes(xmin=order - offset,
                           xmax=order + offset, 
                           ymin=min,
                           ymax=max, fill=sector)) +
    scale_x_continuous(breaks=unique(df$order), labels=unique(df$category))
  
  return(gg)
}

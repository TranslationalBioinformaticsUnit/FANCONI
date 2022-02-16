our.likert=function (items, summary, grouping = NULL, factors = NULL, importance, 
          nlevels = length(levels(items[, 1]))) 
{
  if (!missing(summary)) {
    if (!is.null(grouping) & length(grouping) != nrow(summary)) {
      stop("The length of grouping must be equal to the number of rows in summary.")
    }
    if (!missing(importance)) {
      stop("Pre-summarized data with importance is not currently supported")
    }
    if (is.null(grouping)) {
      r <- list(results = summary, items = NULL, grouping = grouping, 
                nlevels = (ncol(summary) - 1), levels = names(summary[, 
                                                                      2:ncol(summary)]))
    }
    else {
      r <- list(results = cbind(Group = grouping, summary), 
                items = NULL, grouping = grouping, nlevels = (ncol(summary) - 
                                                                1), levels = names(summary[, 2:ncol(summary)]))
    }
    class(r) <- "likert"
    return(r)
  }
  else {
    if (!("data.frame" %in% class(items))) {
      stop(paste0("The items parameter must be a data frame. If trying ", 
                  "to subset a data frame to analyze only one column, try: ", 
                  "items=mydf[,1, drop=FALSE]."))
    }
    if (!all(sapply(items, function(x) "factor" %in% 
                    class(x)))) {
      warning("items parameter contains non-factors. Will convert to factors")
      for (i in 1:ncol(items)) {
        items[, i] <- factor(items[, i], levels = 1:nlevels)
      }
    }
    if (!all(sapply(items, function(x) {
      length(levels(x))
    }) == nlevels)) {
      stop("All items (columns) must have the same number of levels")
    }
    if (!missing(importance)) {
      if (ncol(importance) != ncol(items) | nrow(importance) != 
          nrow(items)) {
        stop("The dimensions of items and importance must be the same")
      }
      if (!all(sapply(importance, function(x) "factor" %in% 
                      class(x)))) {
        warning("importance parameter contains non-factors. Will convert to factors")
        for (i in 1:ncol(importance)) {
          importance[, i] <- factor(importance[, i], 
                                    levels = 1:nlevels)
        }
      }
      if (!all(sapply(importance, function(x) {
        length(levels(x))
      }) == nlevels)) {
        stop("All columns in importance must have the same number of levels")
      }
    }
    if (!is.null(grouping) & !is.null(factors)) {
      stop("Grouping by a grouping column and factor is not supported. \n\t\t\t\t Specify either grouping or factors.")
    }
    lowrange <- 1:ceiling(nlevels/2 - nlevels%%2)
    highrange <- ceiling(nlevels/2 + 1):nlevels
    results <- data.frame()
    if (!is.null(grouping)) {
      if (is.numeric(grouping)) {
        grouping <- as.character(grouping)
      }
      results <- data.frame(Group = rep(unique(grouping), 
                                        each = nlevels), Response = rep(1:nlevels, length(unique(grouping))))
      results2 <- data.frame(Group = rep(unique(grouping), 
                                        each = nlevels), Response = rep(1:nlevels, length(unique(grouping))))
      for (i in 1:ncol(items)) {
        t <- as.data.frame(table(grouping, items[, i]))
        t <- reshape2::dcast(t, Var2 ~ grouping, value.var = "Freq", 
                             add.missing = TRUE)
        t2 <- cbind(Response = t[, 1], apply(t[, 2:ncol(t)], 
                                            2, FUN = function(x) {
                                              x/sum(x) * 100
                                            }))
        t1 <- cbind(Response = t[, 1], apply(t[, 2:ncol(t)], 
                                             2, FUN = function(x) {
                                               x
                                             }))
        t2 <- reshape2::melt(t2)
        t1 <- reshape2::melt(t1)
        results <- merge(results, t2, by.x = c("Group", 
                                              "Response"), by.y = c("Var2", "Var1"),all.x = TRUE)
        names(results)[ncol(results)] <- paste0("Col",i)
        
        results2 <- merge(results2, t1, by.x = c("Group", 
                                              "Response"), by.y = c("Var2", "Var1"),all.x = TRUE)
        names(results2)[ncol(results2)] <- paste0("Col",i)
        
      }
      ##Trequencies
      names(results)[3:ncol(results)] <- names(items)
      results$Response <- factor(results$Response, levels = 1:nlevels, 
                                 labels = levels(items[, i]))
      results <- reshape2::melt(results, id = c("Group", 
                                                "Response"))
      results <- reshape2::dcast(results, Group + variable ~ 
                                   Response)
      results <- as.data.frame(results)
      names(results)[2] <- "Item"
      for (i in 3:ncol(results)) {
        narows <- which(is.na(results[, i]))
        if (length(narows) > 0) {
          results[narows, i] <- 0
        }
      }
      ##Absolut
      names(results2)[3:ncol(results2)] <- names(items)
      results2$Response <- factor(results2$Response, levels = 1:nlevels, 
                                 labels = levels(items[, i]))
      results2 <- reshape2::melt(results2, id = c("Group", 
                                                "Response"))
      results2 <- reshape2::dcast(results2, Group + variable ~ 
                                   Response)
      results2 <- as.data.frame(results2)
      names(results2)[2] <- "Item"
      for (i in 3:ncol(results2)) {
        narows <- which(is.na(results2[, i]))
        if (length(narows) > 0) {
          results2[narows, i] <- 0
        }
      }
    }
    else {
      results <- data.frame(Response = 1:nlevels)
      results2 <- data.frame(Response = 1:nlevels)
      means <- numeric()
      sds <- numeric()
      for (i in 1:ncol(items)) {
        t <- table(items[, i])
        t2 <- (t/sum(t) * 100)
        means[i] <- mean(as.numeric(items[, i]), na.rm = TRUE)
        sds[i] <- sd(as.numeric(items[, i]), na.rm = TRUE)
        results <- cbind(results, as.data.frame(t2)[,2])
        names(results)[ncol(results)] <- names(items)[i]
        results2 <- cbind(results2, as.data.frame(t)[,2])
        names(results2)[ncol(results2)] <- names(items)[i]
      }
      ##freq results
      results <- as.data.frame(t(results))
      names(results) <- levels(items[, 1])
      results <- results[2:nrow(results), ]
      results <- cbind(row.names(results), results)
      names(results)[1] <- "Item"
      row.names(results) <- 1:nrow(results)
      for (i in 2:ncol(results)) {
        narows <- which(is.na(results[, i]))
        if (length(narows) > 0) {
          results[narows, i] <- 0
        }
      }
      ##abs results
      results2 <- as.data.frame(t(results2))
      names(results2) <- levels(items[, 1])
      results2 <- results2[2:nrow(results2), ]
      results2 <- cbind(row.names(results2), results2)
      names(results2)[1] <- "Item"
      row.names(results2) <- 1:nrow(results2)
      for (i in 2:ncol(results2)) {
        narows <- which(is.na(results2[, i]))
        if (length(narows) > 0) {
          results2[narows, i] <- 0
        }
      }
      
      
      if (!is.null(factors)) {
        if (length(factors) != nrow(results)) {
          stop("length(factors) != ncol(items)")
        }
        results <- cbind(results[, 1, drop = FALSE], 
                         Group = factors, results[, 2:ncol(results)])
        results2 <- cbind(results2[, 1, drop = FALSE], 
                         Group = factors, results2[, 2:ncol(results2)])
      }
    }
    r <- list(results = results, results2=results2, items = items, grouping = grouping, 
              factors = factors, nlevels = nlevels, levels = levels(items[, 
                                                                          1]))
    if (!missing(importance)) {
      names(importance) <- names(items)
      r$importance <- likert(importance, grouping = grouping, 
                             nlevels = nlevels)
      class(r) <- c("likert.gap", "likert")
    }
    else {
      class(r) <- "likert"
    }
    return(r)
  }
}


our.summary.likert <- function(object, center=(object$nlevels-1)/2 + 1,
                           ordered=TRUE, ...) {
  if(center < 1.5 | center > (object$nlevels - 0.5) | center %% 0.5 != 0) {
    stop(paste0('Invalid center. Values can range from 1.5 to ', 
                (object$nlevels - 0.5), ' in increments of 0.5'))
  }
  
  if(is.null(object$items)) { # Pre-summarized data
    results <- object$results2
    startCol <- ifelse(is.null(object$grouping), 2, 3)
    if(any(apply(results[,startCol:ncol(results)], 1, FUN=sum) < 1.05)) { 
      # Something a little bigger than 1 to account for rounding errors. 
      # If TRUE, then the percentages range from 0 to 1.
      results[,startCol:ncol(results)] <- 100 * results[,startCol:ncol(results)]
    } #TODO: Could add a check to ensure the sum of all rows equal 1 or 100
    tmp <- t(apply(results[,startCol:ncol(results)], 1, FUN=function(x) {
      tmp <- rep(1:length(x), x)
      return(c(mean=mean(tmp), sd=sd(tmp)))
    }))
    if(center %% 1 == 0) {
      neutral <- results[,(center+(startCol-1))]
    } else {
      neutral <- NA
    }
    if(is.null(object$grouping)) {
      lowCols <- startCol:floor(center + 0.5)
      if(length(lowCols) == 1) {
        low <- results[,lowCols]
      } else {
        low <- apply(results[,lowCols], 1, sum)
      }
      highCols <- (floor(center)+startCol):ncol(results)
      if(length(highCols) == 1) {
        high <- results[,highCols]
      } else {
        high <- apply(results[,highCols], 1, sum)
      }
      results2 <- data.frame(Item=results[,1],
                             low=low,
                             neutral=neutral,
                             high=high,
                             mean=tmp[,1],
                             sd=tmp[,2]  )
    } else {
      lowCols <- startCol:(ceiling(center) + 1)
      if(length(lowCols) == 1) {
        low <- results[,lowCols]
      } else {
        low <- apply(results[,lowCols], 1, sum)
      }
      highCols <- (floor(center)+startCol):ncol(results)
      if(length(highCols) == 1) {
        high <- results[,highCols]
      } else {
        high <- apply(results[,highCols], 1, sum)
      }
      results2 <- data.frame(Group=results[,1],
                             Item=results[,2],
                             low=low,
                             neutral=neutral,
                             high=high,
                             mean=tmp[,1],
                             sd=tmp[,2]  )			
    }
    narows <- which(is.na(results2$low))
    if (length(narows) > 0) {
      results2[narows, ]$low <- 0
    }
    narows <- which(is.na(results2$neutral))
    if (length(narows) > 0) {
      results2[narows, ]$neutral <- 0
    }
    narows <- which(is.na(results2$high))
    if (length(narows) > 0) {
      results2[narows, ]$high <- 0
    }
    return(results2)
  } else {
    results <- object$results2
    items <- object$items
    nlevels <- object$nlevels
    lowrange <- 1 : floor(center - 0.5)
    highrange <- ceiling(center + 0.5) : nlevels
    if(!is.null(object$grouping)) { #Grouping
      results2 = data.frame(Group=rep(unique(results$Group), each=ncol(items)),
                            Item=rep(names(items), length(unique(results$Group))), 
                            low=rep(NA, ncol(items) * length(unique(results$Group))), 
                            neutral=rep(NA, ncol(items) * length(unique(results$Group))),
                            high=rep(NA, ncol(items) * length(unique(results$Group))),
                            mean=rep(NA, ncol(items) * length(unique(results$Group))),
                            sd=rep(NA, ncol(items) * length(unique(results$Group))) )
      for(g in unique(results$Group)) {
        if(length(lowrange) == 1) {
          results2[which(results2$Group == g),]$low <-
            results[results$Group == g, lowrange + 2]
        } else {
          results2[which(results2$Group == g),]$low = apply(
            results[results$Group == g, lowrange + 2], 1, sum)
        }
        if(length(highrange) == 1) {
          results2[which(results2$Group == g),]$high <-
            results[results$Group == g, highrange + 2]
        } else {
          results2[which(results2$Group == g),]$high = apply(
            results[results$Group == g,highrange + 2], 1, sum)
        }
        if(lowrange[length(lowrange)] + 1 != highrange[1]) {
          results2[which(results2$Group == g),]$neutral <- 
            results[results$Group == g, (highrange[1] - 1 + 2)]
        }
        for(i in names(items)) {
          results2[which(results2$Group == g & results2$Item == i), 'mean'] = 
            mean(as.numeric(items[which(object$grouping == g), i]), na.rm=TRUE)
          results2[which(results2$Group == g & results2$Item == i), 'sd'] = 
            sd(as.numeric(items[which(object$grouping == g), i]), na.rm=TRUE)
        }
      }
    } else { #No grouping
      results = data.frame(Response=1:nlevels)
      means = numeric()
      sds = numeric()
      for(i in 1:ncol(items)) {
        t = table(items[,i])
        #t = (t / sum(t) * 100)
        means[i] = mean(as.numeric(items[,i]), na.rm=TRUE)
        sds[i] = sd(as.numeric(items[,i]), na.rm=TRUE)
        results = cbind(results, as.data.frame(t)[2])
        names(results)[ncol(results)] = names(items)[i]
      }
      results = as.data.frame(t(results))
      names(results) = levels(items[,1])
      results = results[2:nrow(results),]
      results2 = data.frame(Item=row.names(results),
                            low=rep(NA, nrow(results)),
                            neutral=rep(NA, nrow(results)),
                            high=rep(NA, nrow(results)),
                            mean=means, 
                            sd=sds)
      if(length(lowrange) == 1) {
        results2$low <- results[,lowrange]
      } else {
        results2$low <- apply(results[,lowrange], 1, sum)
      }
      if(length(highrange) == 1) {
        results2$high <- results[,highrange]
      } else {
        results2$high <- apply(results[,highrange], 1, sum)
      }
      if(lowrange[length(lowrange)] + 1 != highrange[1]) {
        results2$neutral <- results[,(highrange[1] - 1)]
      }
      row.names(results2) = 1:nrow(results2)
      if(!is.null(object$factors)) {
        results2 <- cbind(results2[,1,drop=FALSE],
                          'Factor'=object$factors,
                          results2[,2:ncol(results2),drop=FALSE])
      }
      if(ordered) {
        results2 <- results2[order(results2$high, decreasing=TRUE),]
      }
    }
    
    narows <- which(is.na(results2$low))
    if(length(narows) > 0) {
      results2[narows,]$low <- 0
    }
    narows <- which(is.na(results2$neutral))
    if(length(narows) > 0) {
      results2[narows,]$neutral <- 0
    }
    narows <- which(is.na(results2$high))
    if(length(narows) > 0) {
      results2[narows,]$high <- 0
    }
    
    #Check for factor order in original object and apply to the summary
    if(is.factor(object$results2$Item)) {
      results2$Item <- factor(results2$Item, levels=levels(object$results$Item))
    }
    
    return(results2)
  }
}


##Internal functions needed for plot
label_wrap_mod <- function(value, width = 25) {
  sapply(strwrap(as.character(value), width=width, simplify=FALSE), 
         paste, collapse="\n")
}

abs_formatter <- function(x) {
  return(abs(x))
}

our.likert.bar.plot=function (l, group.order, center = (l$nlevels - 1)/2 + 1, ...) 
{
  opts <- likert.options(...)
  for (i in names(opts)) {
    assign(i, opts[[i]], environment())
  }
  if (center < 1.5 | center > (l$nlevels - 0.5) | center%%0.5 != 
      0) {
    stop(paste0("Invalid center. Values can range from 1.5 to ", 
                (l$nlevels - 0.5), " in increments of 0.5"))
  }
 # ymin <- 0
#  ymax <- 100
  ybuffer <- 5
  lowrange <- 1:floor(center - 0.5)
  highrange <- ceiling(center + 0.5):l$nlevels
  cols <- NULL
  if (!is.null(colors) & length(colors) == l$nlevels) {
    cols <- colors
  }
  else {
    if (!is.null(colors) & length(colors) != l$nlevels) {
      warning("The length of colors must be equal the number of levels.")
    }
    ramp <- colorRamp(c(low.color, neutral.color.ramp))
    ramp <- rgb(ramp(seq(0, 1, length = length(lowrange) + 
                           1)), maxColorValue = 255)
    bamp <- colorRamp(c(neutral.color.ramp, high.color))
    bamp <- rgb(bamp(seq(0, 1, length = length(highrange) + 
                           1)), maxColorValue = 255)
    cols <- NULL
    if (center%%1 != 0) {
      cols <- c(ramp[1:(length(ramp) - 1)], bamp[2:length(bamp)])
    }
    else {
      cols <- c(ramp[1:(length(ramp) - 1)], neutral.color, 
                bamp[2:length(bamp)])
    }
  }
  lsum <- our.summary.likert(l, center = center)
  p <- NULL
  if (!is.null(l$grouping)) {
    lsum$Item <- label_wrap_mod(lsum$Item, width = wrap)
    l$results2$Item <- label_wrap_mod(l$results2$Item, width = wrap)
    lsum$Group <- label_wrap_mod(lsum$Group, width = wrap.grouping)
    results <- l$results2
    results <- reshape2::melt(results, id = c("Group", 
                                              "Item"))
    results$variable <- factor(results$variable, ordered = TRUE)
    if (TRUE | is.null(l$items)) {
      results$Item <- factor(as.character(results$Item), 
                             levels = unique(results$Item), labels = label_wrap_mod(as.character(unique(results$Item)), 
                                                                                    width = wrap), ordered = TRUE)
    }
    else {
      results$Item <- factor(results$Item, levels = label_wrap_mod(names(l$items), 
                                                                   width = wrap), ordered = TRUE)
    }
    ymin <- 0
    if (centered) {
      ymin <- -1*max(lsum$low)
      rows <- which(results$variable %in% names(l$results)[3:(length(lowrange) + 
                                                                2)])
      results[rows, "value"] <- -1 * results[rows, 
                                             "value"]
      if (center%%1 == 0) {
        rows.mid <- which(results$variable %in% names(l$results)[center + 
                                                                   2])
        if (include.center) {
          tmp <- results[rows.mid, ]
          tmp$value <- tmp$value/2 * -1
          results[rows.mid, "value"] <- results[rows.mid, 
                                                "value"]/2
          results <- rbind(results, tmp)
        }
        else {
          results <- results[-rows.mid, ]
        }
      }
      results.low <- results[results$value < 0, ]
      results.high <- results[results$value > 0, ]
      ggplot2.version <- as.integer(unlist(strsplit(as.character(utils::packageVersion("ggplot2")), 
                                                    split = ".", fixed = TRUE)))
      if (ggplot2.version[1] == 2 & ggplot2.version[2] >= 
          2 | ggplot2.version[1] > 2) {
        results.high$variable <- factor(as.character(results.high$variable), 
                                        levels = rev(levels(results.high$variable)))
      }
      p <- ggplot(results, aes(y = value, x = Group, group = variable)) + 
        geom_hline(yintercept = 0) + geom_bar(data = results.low[nrow(results.low):1, 
        ], aes(fill = variable), stat = "identity") + 
        geom_bar(data = results.high, aes(fill = variable), 
                 stat = "identity")
      names(cols) <- levels(results$variable)
      p <- p + scale_fill_manual(legend, breaks = names(cols), 
                                 values = cols, drop = FALSE)
    }
    else {
      ymin <- 0
      p <- ggplot(results, aes(y = value, x = Group, group = variable))
      p <- p + geom_bar(stat = "identity", aes(fill = variable)) + 
        scale_fill_manual(legend, values = cols, breaks = levels(results$variable), 
                          labels = levels(results$variable), drop = FALSE)
    }
    if (plot.percent.low) {
      p <- p + geom_text(data = lsum, y = ymin, aes(x = Group, 
                                                    label = low, group = Item), 
                         size = text.size, hjust = 1, color = text.color)
    }
    if (plot.percent.high) {
      p <- p + geom_text(data = lsum, aes(x = Group, y = max(lsum$high), 
                                          label = high, group = Item), 
                         size = text.size, hjust = -0.2, color = text.color)
    }
    if (plot.percent.neutral & l$nlevels%%2 == 1 & include.center) {
      if (centered) {
        p <- p + geom_text(data = lsum, y = 0, aes(x = Group, 
                                                   group = Item, label = paste0(round(neutral), 
                                                                                "%")), size = text.size, hjust = 0.5, 
                           color = text.color)
      }
      else {
        lsum$y <- lsum$low + (lsum$neutral/2)
        p <- p + geom_text(data = lsum, aes(x = Group, 
                                            y = y, group = Item, label = paste0(round(neutral), 
                                                                                "%")), size = text.size, hjust = 0.5, 
                           color = text.color)
      }
    }
    if (FALSE & plot.percents) {
      warning("plot.percents is not currenlty supported for grouped analysis.")
    }
    p <- p + coord_flip() + ylab("Absolut frequencies") + xlab("") + 
      theme(axis.ticks = element_blank(), strip.background = element_rect(fill = panel.strip.color, 
                                                                          color = panel.strip.color))
    if (is.null(panel.arrange)) {
      p <- p + facet_wrap(~Item)
    }
    else if (panel.arrange == "v") {
      p <- p + facet_wrap(~Item, ncol = 1)
    }
    else if (panel.arrange == "h") {
      p <- p + facet_wrap(~Item, nrow = 1)
    }
    if (!missing(group.order)) {
      p <- p + scale_x_discrete(limits = rev(group.order), 
                                drop = FALSE)
    }
  }
  else {
    factor.mapping <- NULL
    if (!is.null(l$factors)) {
      factor.mapping <- l$results[, 1:2]
      names(factor.mapping)[2] <- "Factor"
      results <- reshape2::melt(l$results2[, -2], id.vars = "Item")
    }
    else {
      results <- reshape2::melt(l$results2, id.vars = "Item")
    }
    if (ordered & is.null(results$factor)) {
      #order <- lsum[order(lsum$high), "Item"]
      order <- lsum[order(lsum$Item, decreasing=TRUE), "Item"]
      results$Item <- factor(results$Item, levels = order)
    }
    ymin <- 0
    if (centered) {
      ymin <- -1*max(lsum$low)
      rows <- which(results$variable %in% names(l$results2)[2:(length(lowrange) + 
                                                                1)])
      results[rows, "value"] <- -1 * results[rows, 
                                             "value"]
      if (center%%1 == 0) {
        rows.mid <- which(results$variable %in% names(l$results)[center + 
                                                                   1])
        if (include.center) {
          tmp <- results[rows.mid, ]
          tmp$value <- tmp$value/2 * -1
          results[rows.mid, "value"] <- results[rows.mid, 
                                                "value"]/2
          results <- rbind(results, tmp)
        }
        else {
          results <- results[-rows.mid, ]
        }
      }
      if (!is.null(factor.mapping)) {
        results$order <- 1:nrow(results)
        results <- merge(results, factor.mapping, by = "Item", 
                         all.x = TRUE)
        results <- results[order(results$order), ]
        results$order <- NULL
      }
      results.low <- results[results$value < 0, ]
      results.high <- results[results$value > 0, ]
      p <- ggplot(results, aes(y = value, x = Item, group = Item)) + 
        geom_hline(yintercept = 0) + geom_bar(data = results.low[nrow(results.low):1, 
        ], aes(fill = variable), stat = "identity") + 
        geom_bar(data = results.high, aes(fill = variable), 
                 stat = "identity")
      names(cols) <- levels(results$variable)
      p <- p + scale_fill_manual(legend, breaks = names(cols), 
                                 values = cols, drop = FALSE)
    }
    else {
      if (!is.null(factor.mapping)) {
        results$order <- 1:nrow(results)
        results <- merge(results, factor.mapping, by = "Item", 
                         all.x = TRUE)
        results <- results[order(results$order), ]
        results$order <- NULL
      }
      p <- ggplot(results, aes(y = value, x = Item, group = Item))
      p <- p + geom_bar(stat = "identity", aes(fill = variable))
      p <- p + scale_fill_manual(legend, values = cols, 
                                 breaks = levels(results$variable), labels = levels(results$variable), 
                                 drop = FALSE)
    }
    if (plot.percent.low) {
      p <- p + geom_text(data = lsum, y = ymin, aes(x = Item, 
                                                    label = low), size = text.size, 
                         hjust = 1, color = text.color)
    }
    if (plot.percent.high) {
      p <- p + geom_text(data = lsum, y = max(lsum$high), aes(x = Item, 
                                                   label = high), size = text.size, 
                         hjust = -0.2, color = text.color)
    }
    if (plot.percent.neutral & l$nlevels%%2 == 1 & include.center & 
        !plot.percents) {
      if (centered) {
        p <- p + geom_text(data = lsum, y = 0, aes(x = Item, 
                                                   label = paste0(round(neutral), "%")), 
                           size = text.size, hjust = 0.5, color = text.color)
      }
      else {
        lsum$y <- lsum$low + (lsum$neutral/2)
        p <- p + geom_text(data = lsum, aes(x = Item, 
                                            y = y, label = paste0(round(neutral), "%")), 
                           size = text.size, hjust = 0.5, color = text.color)
      }
    }
    if (plot.percents) {
      center.label <- ""
      if (center%%1 == 0) {
        center.label <- names(l$results)[center + 1]
      }
      lpercentpos <- ddply(results[results$value > 0, ], 
                           .(Item), transform, pos = cumsum(value) - 0.5 * 
                             value)
      p <- p + geom_text(data = lpercentpos[lpercentpos$variable != 
                                              center.label, ], aes(x = Item, y = pos, label = paste0(round(value), 
                                                                                                     "%")), size = text.size, color = text.color)
      lpercentneg <- results[results$value < 0, ]
      if (nrow(lpercentneg) > 0) {
        lpercentneg <- lpercentneg[nrow(lpercentneg):1, 
        ]
        lpercentneg$value <- abs(lpercentneg$value)
        lpercentneg <- ddply(lpercentneg, .(Item), transform, 
                             pos = cumsum(value) - 0.5 * value)
        lpercentneg$pos <- lpercentneg$pos * -1
        p <- p + geom_text(data = lpercentneg[lpercentneg$variable != 
                                                center.label, ], aes(x = Item, y = pos, label = paste0(round(abs(value)), 
                                                                                                       "%")), size = text.size, color = text.color)
      }
      lpercentneutral <- results[results$variable == center.label, 
      ]
      if (nrow(lpercentneutral) > 0) {
        p <- p + geom_text(data = lpercentneutral, aes(x = Item, 
                                                       y = 0, label = paste0(round(abs(value * 2)), 
                                                                             "%")), size = text.size, color = text.color)
      }
    }
    p <- p + coord_flip() + ylab("Absolut frequency") + xlab("") + 
      theme(axis.ticks = element_blank())
    if (!missing(group.order)) {
      p <- p + scale_x_discrete(limits = rev(group.order), 
                                labels = label_wrap_mod(rev(group.order), width = wrap), 
                                drop = FALSE)
    }
    else {
      p <- p + scale_x_discrete(breaks = l$results2$Item, 
                                labels = label_wrap_mod(l$results2$Item, width = wrap), 
                                drop = FALSE)
    }
  }
  p <- p + scale_y_continuous(labels = abs_formatter, limits = c(ymin - 
                                                                   ybuffer, max(lsum$high) + ybuffer))
  p <- p + theme(legend.position = legend.position)
  attr(p, "item.order") <- levels(results$Item)
  class(p) <- c("likert.bar.plot", class(p))
  return(p)
}
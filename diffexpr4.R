library( reshape2 );
library( ggplot2 );
library( scales );

require( edgeR );
require( DESeq2 );
#library( dplyr );
#library( tidyr );

# TODO:
# - getLengthNormCounts
# - boxplot of FC along AvgExpr to see shift from zero?
# - plot number of DE genes vs number of genes per chromosome

getRawCountsExp <- function( x ) UseMethod( "getRawCountsExp" );
getNormalizedExp <- function( x ) UseMethod( "getNormalizedExp" );
getRawCounts <- function( x ) UseMethod( "getRawCounts" );
getCounts <- function( x ) UseMethod( "getCounts" );
getTags <- function( x ) UseMethod( "getTags" );
getDesign <- function( x ) UseMethod( "getDesign" );
getContrasts <- function( x ) UseMethod( "getContrasts" );
getGroups <- function( x ) UseMethod( "getGroups" );
getTagsAnnotations <- function( x ) UseMethod( "getTagsAnnotations" );
getSampleSheet <- function( x ) UseMethod( "getSampleSheet" );
getFormula <- function( x ) UseMethod( "getFormula" );
dendrogramPlot <- function( x, ... ) UseMethod( "dendrogramPlot" );
getEdgeR <- function( x ) UseMethod( "getEdgeR" );
edgeRFit <- function( x, ... ) UseMethod( "edgeRFit" );
getVoom <- function( x ) UseMethod( "getVoom" );
eBayesFit <- function( x ) UseMethod( "eBayesFit" );
gsaKEGG <- function( x, ... ) UseMethod( "gsaKEGG" );
getDESeq2 <- function( x ) UseMethod( "getDESeq2" );
nBinomWaldFit <- function( x ) UseMethod( "nBinomWaldFit" );
deTable <- function( x, ... ) UseMethod( "deTable" );
writeReport <- function( x, ... ) UseMethod( "writeReport" );
plotDiag <- function( x, ... ) UseMethod( "plotDiag" );
pairedCountsPlot <- function( x, ... ) UseMethod( "pairedCountsPlot" );
countsDensityPlot <- function( x, ... ) UseMethod( "countsDensityPlot" );
getSubjectVar <- function( x ) UseMethod( "getSubjectVar" );
getGroupVar <- function( x ) UseMethod( "getGroupVar" );
getSampleVar <- function( x ) UseMethod( "getSampleVar" );
getTagVar <- function( x ) UseMethod( "getTagVar" );
getCountVar <- function( x ) UseMethod( "getCountVar" );
getAlignmentSummary <- function( x ) UseMethod( "getAlignmentSummary" );
countsData <- function( x, ... ) UseMethod( "countsData" );
getAvgCountVar <- function( x ) UseMethod( "getAvgCountVar" );
getDiffCountVar <- function( x ) UseMethod( "getDiffCountVar" );
deData <- function( x, ... ) UseMethod( "deData" );
selectFeatures <- function( x, ... ) UseMethod( "selectFeatures" );

diffTagsPlot <- function( x, ... ) {
  warning( "diffTagsPlot: obsoleted; rename to contrastsDECountPlot" )
  contrastsDECountPlot( x, ... );
}

boxplotCounts <- function( x, ... ) {
  warning( "boxplotCounts: obsoleted; rename to countsBoxPlot" );
  countsBoxPlot( x, ... );
}

plotPairedCounts <- function( x, ... ) {
  warning( "plotPairedCounts: obsoleted; rename to pairedCountsPlot." );
  pairedCountsPlot( x, ... );
}

plotTopTags <- function( x, ... )  {
  warning( "plotTopTags: obsoleted; rename to topTagsPlot." );
  topTagsPlot( x, ... );
}

plot.DETable <- function( x, volcano = FALSE, ... ) {
  warning( "plot.DETable: obsoleted; rename to diffExpPlot or volcanoPlot." );
  if( volcano == TRUE ) {
    volcanoPlot( x, ... );
  } else {
    diffExpPlot( x, ... );
  }
}

plot.CountsExp <- function( x, ... ) {
  warning( "plot.CountsExp: obsoleted; rename to countsDensityPlot." );
  countsDensityPlot( x, ... );
}

plot.RawCountsExp <- function( x, ... ) {
  warning( "plot.CountsExp: obsoleted; rename to countsDensityPlot." );
  countsDensityPlot( x, ... );
}

## -----------------------------------------------------------------------

tableWriter <- function( x, names, formats = c( "tab.gz", "xls" ), ... ) {
  format2writer <- list(
    "-" = function( x, names, sep = "\t", row.names = FALSE, ... ) {
      cat( n, "\n" );
      write.table( head( x ), row.names = row.names, sep = sep, ... );
    },
    "tab" = function( x, names, sep = "\t", row.names = FALSE, ... ) {
      n <- paste0( paste( names, collapse = "/" ), ".tab" );
      dir.create( path = dirname( n ), recursive = TRUE, showWarnings = FALSE );
      write.table( x, file = n, row.names = row.names, sep = sep, ... );
    },
    "tab.gz" = function( x, names, sep = "\t", row.names = FALSE, ... ) {
      n <- paste0( paste( names, collapse = "/" ), ".tab.gz" );
      dir.create( path = dirname( n ), recursive = TRUE, showWarnings = FALSE );
      write.table( x, file = gzfile( n ), row.names = row.names, sep = sep, ... );
    },
    "xlsx" = function( x, names, row.names = TRUE, ... ) {
      library( xlsx );
      n <- paste0( paste( names, collapse = "/" ), ".xlsx" );
      dir.create( path = dirname( n ), recursive = TRUE, showWarnings = FALSE );
      write.xlsx( x = x, file = n, sheetName = paste( names, collapse = "/" ), row.names = row.names );
    },
    "xls" = function( x, names, ... ) {
      library( WriteXLS );
      n <- paste0( paste( names, collapse = "/" ), ".xls" );
      dir.create( path = dirname( n ), recursive = TRUE, showWarnings = FALSE );
      WriteXLS( x = c( "x" ), ExcelFileName = n, ... );
    }
  )

  for( f in formats ) {
    tryCatch(
      format2writer[[ f ]]( x, names, ... ),
      error = function( e ) {
        warning( paste0( "tableWriter (format='", f, "') failed: ", e ) )
      }
    );
  }
}

writeReport.data.frame <- function( x, names = c(), writer = tableWriter, ... ) {
  writer( x, names, ... )
}

writeReport.matrix <- function( x, names = c(), writer = tableWriter, ... ) {
  writer( as.data.frame( x ), names, ... )
}

## -----------------------------------------------------------------------
## internal utils
## -----------------------------------------------------------------------

.getGroupName <- function( x, index ) {
  groups <- getGroups( x );
  if( is.null( groups ) ) {
    NULL
  } else {
    if( !is.null( getFormula( x ) ) ) {
      nms <- labels( terms( getFormula( x ) ) );
      nms <- nms[ !is.na( match( nms, names( groups ) ) ) ];
    } else {
      nms <- names( groups );
      nms <- nms[ sapply( nms, function( nm ) length( unique( getSampleSheet( x )[[ nm ]] ) ) ) > 1 ];
    }
    if( index > length( nms ) ) {
      NULL
    } else {
      nms[[ index ]];
    }
  }
}

## -----------------------------------------------------------------------
## CountsExp
## -----------------------------------------------------------------------

getRawCounts.CountsExp <- function( x ) getRawCounts( getRawCountsExp( x ) );
getGroups.CountsExp <- function( x ) getGroups( getRawCountsExp( x ) );
getSampleSheet.CountsExp <- function( x ) getSampleSheet( getRawCountsExp( x ) );
getFormula.CountsExp <- function( x ) getFormula( getRawCountsExp( x ) );
getDesign.CountsExp <- function( x ) getDesign( getRawCountsExp( x ) );
getContrasts.CountsExp <- function( x ) getContrasts( getRawCountsExp( x ) );
getTagsAnnotations.CountsExp <- function( x ) getTagsAnnotations( getRawCountsExp( x ) );
getTags.CountsExp <- function( x ) rownames( getCounts( x ) );
getSubjectVar.CountsExp <- function( x ) getSubjectVar( getRawCountsExp( x ) );
getGroupVar.CountsExp <- function( x ) getGroupVar( getRawCountsExp( x ) );
getSampleVar.CountsExp <- function( x ) getSampleVar( getRawCountsExp( x ) );
getTagVar.CountsExp <- function( x ) getTagVar( getRawCountsExp( x ) );

countsData.CountsExp <- function( x, addTagsAnnotations = FALSE, geneLenVar = "meanSumExonsLengths", ... ) {
  tagVar <- getTagVar( x );
  sampleVar <- getSampleVar( x );
  countVar <- getCountVar( x );

  if( countVar == "log2CPM" && !is.na( geneLenVar ) ) {
    addTagsAnnotations <- TRUE;
  }

  d <- countsData( getRawCountsExp( x ), addTagsAnnotations = addTagsAnnotations, ... );

  dd <- getCounts( x ) %>% as.data.frame %>% add_rownames( var = tagVar );
  dd <- dd %>% gather_( sampleVar, countVar, setdiff( colnames( dd ), tagVar ) );
  dd[[ tagVar ]] <- factor( dd[[ tagVar ]], levels = levels( d[[ tagVar ]] ) );
  dd[[ sampleVar ]] <- factor( dd[[ sampleVar ]], levels = levels( d[[ sampleVar ]] ) );

  d <- left_join( d, dd, by = c( tagVar, sampleVar ) );

  if( countVar == "log2CPM" && !is.na( geneLenVar ) && ( geneLenVar %in% colnames( d ) ) ) {
    d[[ "log2RPKM" ]] <- d[[ countVar ]] - log2( d[[ geneLenVar ]] / 1000 );
  }

  d
}

writeReport.CountsExp <- function( x, names = c(), ... ) {
  tagVar <- getTagVar( x );
  sampleVar <- getSampleVar( x );

  f <- function( countVar, ... ) {
    d <- countsData( x ) %>%
      select_( .dots = c( tagVar, sampleVar, countVar ) ) %>%
      spread_( sampleVar, countVar );

    writeReport( d, names = c( names, report = countVar ), ... );
  }

  f( getCountVar( x ), ... );
}

# --- plots should go below ---

countsDensityPlot.CountsExp <- function( x, ..., countScale = "continuous" ) {
  countsDensityPlot.general( x, ..., countScale = countScale );
}

pairedCountsPlot.CountsExp <- function( x, ..., countScale = "continuous" ) {
  pairedCountsPlot.general( x, ..., countScale = countScale );
}

plotMDS.CountsExp <- function( x,
                               colorGroup = .getGroupName( x, 1 ),
                               ...
) {
  counts <- getCounts( x );
  groups <- getGroups( x );
  colors <- groups[[ colorGroup ]]$sf[ colnames( counts ), ]$color;

  plotMDS( counts, col = colors, ... );
  grid();
}

rawPcaData <- function( x, dimNum = 5, countVar = getCountVar( x ), sampleVar = getSampleVar( x ), tagVar = getTagVar( x ), ... ) {
  library( FactoMineR );

  d <- countsData( x, ... ) %>%
    dplyr::select_( .dots = c( tagVar, sampleVar, countVar ) ) %>%
    spread_( tagVar, countVar ) %>%
    as.data.frame();
  rownames( d ) <- d[[ sampleVar ]];
  d <- d[ setdiff( colnames( d ), sampleVar ) ];
  pca <- PCA( d, graph = FALSE, ncp = dimNum );

  list(
    pca = pca
  )
}

pcaDimsData <- function( x, dims = c( 1, 2 ), frame = getSampleSheet( x ), sampleVar = getSampleVar( x ), ... ) {
  pd <- rawPcaData( x, ... );
  d <- pd$pca$ind$coord[ , dims, drop = FALSE ];

  dd <- do.call( rbind, lapply( 1:( ncol( d ) - 1 ), function( row ) {
    do.call( rbind, lapply( ( row + 1 ):ncol( d ), function( col ) {
      data.frame(
        row = colnames( d )[[row]],
        col = colnames( d )[[col]],
        rowD = d[ , row ],
        colD = d[ , col ],
        .id = rownames( d )
      );
    } ) )
  } ) );
  dd <- merge( dd, frame, by.x = ".id", by.y = "row.names" ) %>% dplyr::select( -.id );
  dd
}

pcaDimsPlot <- function( x, labelVar = getSampleVar( x ),
  colorVar = .getGroupName( x, 1 ), shapeVar = NULL, lineVar = NULL,
  pointSize = 3, labelSize = 3, pointAlpha = 0.5,
  ...
) {
  dd <- pcaDimsData( x, ... );

  if( !is.null( shapeVar ) ) dd[[ shapeVar ]] <- factor( dd[[ shapeVar ]] );
  if( !is.null( colorVar ) ) dd[[ colorVar ]] <- factor( dd[[ colorVar ]] );

  p <- ggplot( dd, aes_string(
    x = "colD",
    y = "rowD",
    label = labelVar,
    color = colorVar,
    shape = shapeVar
  ) );
  if( !is.null( lineVar ) ) {
    p <- p + geom_line( aes_string( group = lineVar ), linetype = 2, color = "gray", alpha = 0.3 );
  }
  if( is.null( labelVar ) ) {
    p <- p + geom_point( size = pointSize, alpha = pointAlpha );
  } else {
    p <- p + geom_text( size = labelSize );
  }
  if( !is.null( shapeVar ) ) {
    p <- p + scale_shape_manual( values = 1:nlevels( dd[[ shapeVar ]] ) %% 6 );
  }
  p <- p + facet_grid( row ~ col, scale = "free" ) + xlab( "" ) + ylab( "" );
  p <- p + theme_bw();
  p
}

pcaData.CountsExp <- function( x, frame = getSampleSheet( x ), groupNames = NULL, maxLevelsNum = nrow( frame ) - 1, ... ) {
  pd <- rawPcaData( x, ... );
  pca <- pd$pca;

  if( is.null( groupNames ) ) groupNames <- colnames( frame );
  groupNames <- groupNames[ sapply( groupNames, function( gn ) class( frame[[ gn ]] ) ) %in% c( "factor", "character" ) ];
  groupNames <- groupNames[ sapply( groupNames, function( gn ) nlevels( factor( frame[[ gn ]] ) ) ) < maxLevelsNum ];
  groupNames <- groupNames[ sapply( groupNames, function( gn ) nlevels( factor( frame[[ gn ]] ) ) ) > 1 ];

  frame <- frame[ groupNames ];
  groupValues <- unlist( lapply( groupNames, function( gn ) {
    as.character( if( is.factor( frame[ gn ] ) ) levels( frame[[ gn ]] ) else unique( frame[[ gn ]] ) )
  } ) );
  frame <- data.frame( apply( frame, c(1,2), as.character ), stringsAsFactors = FALSE );
  frame[[ getSampleVar( x ) ]] <- rownames( frame );
  d <- melt( pca$ind$coord );
  colnames( d ) <- c( getSampleVar( x ), "Dim", "DimValue" );
  d <- merge( d, frame, by = getSampleVar( x ) );
  d <- melt( d, measure.vars = groupNames, value.name = "GroupValue", variable.name = "Group" );
  d$GroupValue <- factor( d$GroupValue, levels = groupValues );

  l <- list(
    d = d,
    pca = pca
  );
  class( l ) <- c( "PcaData", class( l ) );
  l
}

pcaEigenVarData <- function( x, ... ) {
  pd <- rawPcaData( x, ... );
  d <- pd$pca$eig;
  d <- data.frame( eigenvalue = d$eigenvalue, percentage = d[[ "percentage of variance" ]], cumulative = d[[ "cumulative percentage of variance" ]] );
  d$dim <- 1:nrow( d );
  d
}

pcaEigenVarPlot <- function( x, ... ) {
  d <- pcaEigenVarData( x, ... );
  #d <- melt( d, id.vars = c( "dim", "eigenvalue" ), variable.name = "varianceType", value.name = "value" );
  #ggplot( d, aes( x = dim, y = value, color = varianceType, label = dim ) ) + geom_line() + geom_text() + theme_bw()
  ggplot( d, aes( x = cumulative, y = percentage, label = dim ) ) +
    geom_line( color = "gray", linetype = 2 ) +
    geom_text( angle = -45 ) +
    theme_bw();
}

pcaPlot <- function( x, labelVar = TRUE, labelSize = 3, position = "jitter", ... ) {
  if( is.logical( labelVar ) && labelVar ) labelVar <- getSampleVar( x );
  l <- pcaData.CountsExp( x, ... );
  d <- l$d;

  p <- ggplot( d, aes_string(
    x = "GroupValue",
    y = "DimValue",
    label = labelVar
  ) );
  if( is.null( labelVar ) ) {
    p <- p + geom_point( position = position );
  } else {
    p <- p + geom_text( size = labelSize );
  }
  p <- p +
    facet_grid( Dim ~ Group, scale = "free" ) +
    xlab( "Groups and Group Values") +
    ylab( "Principal Components Values" ) +
    theme_bw() +
    theme( axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1 ) );
  p
}

dendrogramPlot.CountsExp <- function( x,
                                      colorGroup = .getGroupName( x, 1 ),
                                      colorGroups = c( colorGroup ),
                                      labelVar = NULL,
                                      ...
) {
  counts <- getCounts( x );
  groups <- getGroups( x );

  if( !is.null( labelVar ) ) {
    lbls <- groups[[ labelVar ]]$sf[ colnames( counts ), ]$label;
    colnames( counts ) <- lbls;
  }

  d <- dist( t( counts ) );
  hc <- hclust( d );
  dgr <- as.dendrogram( hc );

  if( 0 ) {
    colors <- groups[[ colorGroup ]]$sf[ colnames( counts ), ]$color;

    library( ggdendro );

    #plot( hclust( d ), cex = 0.6 )

    #ggdendrogram( hc, rotate = TRUE );

    dgrData <- dendro_data( dgr );
    labs <- label( dgrData );
    labs$color <- groups[[ colorGroup ]]$sf[ as.character( labs$label ), ]$label;

    ggplot( segment( dgrData ) ) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text( data = labs, aes( label = label, x = x, y = 0, color = color ), hjust = 1 ) +
      coord_flip() +
      scale_y_continuous(expand = c(0.3, 0)) +
      theme_bw();
  } else {
    library( dendextend );

    nms <- rev( colorGroups );
    clrs <- do.call( cbind, lapply( nms, function( g ) groups[[ g ]]$sf[ order.dendrogram( dgr ), ]$color ) );
    colnames( clrs ) <- nms;

    dgr %>% hang.dendrogram( hang_height = 10 ) %>% plot( ... )
    colored_bars( colors = clrs, dend = dgr )
  }
}

groupsSummary <- function( countsExp ) {
  counts <- getCounts( countsExp );
  groups <- getGroups( countsExp );
  l <- lapply( names( groups ), function( groupName ) {
    d <- data.frame( groups[[ groupName ]]$sf$label );
    colnames( d )[[1]] <- groupName;
    d
  } );
  d <- do.call( cbind, l );
  rownames( d ) <- colnames( counts );
  #  names( dimnames( d ) ) <- c( names( dimnames( counts ) )[[2]], "group" );
  d
}

if( 0 ) { # switched off on 20160605
  plotSplitMDS.CountsExp <- function( countsExp,
                                      groups = getGroups( countsExp ),
                                      colorGroup = if( is.null( groups ) ) NULL else names( groups )[[1]],
                                      tagsGroups = NULL,
                                      top = 500,
                                      minTags = 30,
                                      dim.plot = c( 1, 2 ), ndim = max( dim.plot ),
                                      gene.selection = "pairwise",
                                      showLabels = FALSE,
                                      ...
  ) {
    x <- getCounts( countsExp );

    nsamples <- ncol( x );
    stopifnot( nsamples >= 3 );
    stopifnot( ndim >= 2 );
    stopifnot( nsamples >= ndim );

    cn <- colnames( x );
    bad <- rowSums(is.finite(x)) < nsamples;
    if( any( bad ) ) x <- x[ !bad, , drop = FALSE ];

    calcMDS <- function( x, top ) {
      nprobes <- nrow( x );
      if( nprobes >= ndim && nprobes >= minTags ) {
        top <- min( top, nprobes );

        dd <- matrix( 0, nrow = nsamples, ncol = nsamples, dimnames = list(cn, cn) );
        topindex <- nprobes - top + 2;
        if (gene.selection == "pairwise") {
          for (i in 2:(nsamples))
            for (j in 1:(i - 1))
              dd[i, j] = sqrt(mean(sort.int((x[, i] - x[, j])^2, partial = topindex)[topindex:nprobes]))
        } else {
          s <- rowMeans((x - rowMeans(x))^2)
          q <- quantile(s, p = (topindex - 1.5)/(nprobes - 1))
          x <- x[s >= q, ]
          for (i in 2:(nsamples))
            dd[i, 1:(i - 1)] = sqrt(colMeans((x[, i] - x[, 1:(i - 1), drop = FALSE])^2))
        }
        a1 <- cmdscale(as.dist(dd), k = ndim)
        d <- data.frame( x = a1[, dim.plot[1]], y =- a1[, dim.plot[2]] );
        d$sample <- rownames( d );
        d
      } else {
        data.frame();
      }
    }

    if( is.null( tagsGroups ) ) {
      d <- calcMDS( x, top );
      d$group <- "all";
    } else {
      l <- lapply( split( 1:nrow( x ), tagsGroups ), function( rows ) {
        calcMDS( x[ rows, , drop = FALSE ], top );
      } );
      d <- do.call( rbind, lapply( names( l ), function( n ) { d <- l[[n]]; if( nrow( d ) > 0 ) d$group <- n; d } ) );
    }

    d$color <- groups[[ colorGroup ]]$sf[ as.character( d$sample ), ]$label;

    p <- ggplot( d, aes( x = x, y = y, label = sample, color = color, shape = color ) ) + facet_wrap( ~ group ) + theme_bw();
    if( showLabels ) {
      p <- p + geom_text( size = 3 );
    } else {
      p <- p + geom_point( size = 2, alpha = 0.7 );
    }
    p
  }
}

#fromVar toVar fromValue toValue
#sampleId color s1 red
mergeProps <- function( x, d, p, .propMap = NULL ) {
  lJoin <- c();
  if( is.null( p$fromCol ) ) {
    lJoin[[ p$fromVar ]] <- p$fromVar;
  } else {
    lJoin[[ p$fromCol ]] <- p$fromVar;
  }
  lSelect <- c( p$fromVar );
  if( is.null( p$toCol ) ) {
    lSelect[[ p$toVar ]] <- p$toVar;
  } else {
    lSelect[[ p$toCol ]] <- p$toVar;
  }
  pm <- .propMap %>% dplyr::select_( .dots = lSelect );
  d %>% inner_join( pm, by = lJoin );
}

### unittest {
.d <- data.frame( SampleId = c( "a", "b", "c", "b" ) );
.map <- data.frame( sampleId = c( "b", "a", "c" ), color = c( "r", "g", "b" ), shape = 1:3 );
stopifnot( identical(
  mergeProps( x = NULL, d = .d, list( fromCol = "SampleId", fromVar = "sampleId", toVar = "color", toCol = "Color" ), .propMap = .map ),
  data.frame( SampleId = c( "a", "b", "c", "b" ), Color = c( "g", "r", "b", "r" ) )
) );
### } # unittest

countsDensityData <- function( x,
  tagsAnnotations = NULL,
  zeroCount = 0,
  filterFun = function( d ) d,
  ...
) {
  countVar <- getCountVar( x );

  ### ---- obsleting (start)
  if( 0 ) { # obsoleting in 20160504
    sampleVar <- getSampleVar( x );
    tagVar <- getTagVar( x );

    counts <- getCounts( x );
    d <- melt( counts, value.name = countVar );
    stopifnot( sampleVar %in% colnames( d ) );
    stopifnot( tagVar %in% colnames( d ) );

    if( is.numeric( d[[ sampleVar ]] ) ) d[[ sampleVar ]] <- factor( d[[ sampleVar ]] );
    if( is.numeric( d[[ tagVar ]] ) ) d[[ tagVar ]] <- factor( d[[ tagVar ]] );

    ### ----- add group annotations -----
    groups <- getGroups( x );
    for( gn in names( groups ) ) {
      if( !is.null( gn ) && gn != sampleVar ) {
        if( 1 ) {
          pm <- groups[[ gn ]]$sf;
          pm[[ sampleVar ]] <- rownames( pm );
          d <- mergeProps( x = x, d = d, list( fromVar = sampleVar, toVar = "level", toCol = gn ), pm );
        } else {
          d[[ gn ]] <- groups[[ gn ]]$sf$level[ d[[ sampleVar ]] ];
        }
      }
    }

    ### ----- add tag annotations -----
    if( is.logical( tagsAnnotations ) && tagsAnnotations == TRUE ) {
      tagsAnnotations <- getTagsAnnotations( x );
      if( "chromosome_name" %in% colnames( tagsAnnotations ) ) {
        dd <- as.data.frame( table( chromosome_name = tagsAnnotations$chromosome_name, useNA = "always" ) );
        dd$chromosome_name_count <- with( dd, paste0( chromosome_name, " [", Freq, "]" ) );
        dd$chromosome_name_count <- factor( dd$chromosome_name_count, levels = dd$chromosome_name_count );
        tagsAnnotations$chromosome_name_count <- dd$chromosome_name_count[ match( tagsAnnotations$chromosome_name, dd$chromosome_name ) ];
      }
    }

    if( !is.null( tagsAnnotations ) ) {
      d <- data.frame( d, tagsAnnotations[ d[[ tagVar ]], ] );
    }
  }
  ### ---- obsleting (end)

  d <- countsData( x, addSampleSheet = TRUE, addTagsAnnotations = TRUE, ... );
  d[[ countVar ]][ d[[ countVar ]] == 0 ] <- zeroCount;

  filterFun( d );
}

countsDensityPlot.general <- function( x,
  colorGroup = .getGroupName( x, 1 ), # obsoleting on 20161101
  colorVar = colorGroup,
  linetypeGroup = .getGroupName( x, 2 ),
  linetypeVar = linetypeGroup, # obsoleting on 20161101
  countVar = getCountVar( x ),
  countScale = "log10",
  ...
) {
  d <- countsDensityData( x, ... );

  groups <- getGroups( x );

  if( countScale == "log10" ) {
    d <- subset( d, d[[ countVar ]] > 0 );
  }

  if( nrow( d ) > 0 ) {
    p <- ggplot( d, aes_string(
      x = countVar,
      color = if( !is.null( colorVar ) ) colorVar else NULL,
      linetype = if( !is.null( linetypeVar ) ) linetypeVar else NULL,
      group = getSampleVar( x )
    ) ) +
      theme_bw();
    p <- p + stat_density( geom = "path", position = "identity", kernel = "rectangular", alpha = 0.5 );
    if( !is.null( colorVar ) ) {
      colorValues <- groups[[ colorVar ]]$df$color;
      names( colorValues ) <- groups[[ colorVar ]]$df$level;
      p <- p + scale_color_manual( name = groups[[ colorVar ]]$name, values = colorValues );
    }
    if( !is.null( linetypeVar ) ) {
      linetypeValues <- as.factor( groups[[ linetypeVar ]]$df$linetype );
      names( linetypeValues ) <- groups[[ linetypeVar ]]$df$level;
      p <- p + scale_linetype_manual( name = groups[[ linetypeVar ]]$name, values = linetypeValues );
    }
    if( countScale == "log10" ) {
      p <- p + scale_x_log10();
    }
    p
  } else {
    NULL
  }
}

pairedCountsData <- function( x,
  groupLevels = as.character( head( ( getGroups( x ) )[[ getGroupVar( x ) ]]$df$level, 2 ) ),
  countVar = getCountVar( x ),
  ...
) {
  subjectVar <- getSubjectVar( x );
  groupVar <- getGroupVar( x );

  #df <- countsDensityData( x, ... );
  df <- countsData( x, addSampleSheet = TRUE, addTagsAnnotations = TRUE, addDePattern = TRUE, ... );

  dx <- df[ df[[ groupVar ]] %in% groupLevels[[1]], ];
  dy <- df[ df[[ groupVar ]] %in% groupLevels[[2]], ];
  rn <- range( dx[[ countVar ]], dy[[ countVar ]], na.rm = TRUE, finite = TRUE );
  df <- merge( dx, dy, by = c( subjectVar, getTagVar( x ) ) );

  df
}

pairedCountsPlot.general <- function( x, ...,
  groupLevels = as.character( head( ( getGroups( x ) )[[ getGroupVar( x ) ]]$df$level, 2 ) ),
  countVar = getCountVar( x ),
  colorGroup = NULL, # obsoleting, will be removed 20160901
  colorVar = colorGroup,
  countScale = "log10",
  showFit = TRUE
) {
  df <- pairedCountsData( x, groupLevels = groupLevels, countVar = countVar, ... );

  subjectVar <- getSubjectVar( x );
  groupVar <- getGroupVar( x );
  groups <- getGroups( x );

  xCountVar <- paste0( countVar, ".x" );
  yCountVar <- paste0( countVar, ".y" );

  rows <- !is.na( df[[ xCountVar ]] ) & !is.na( df[[ yCountVar ]] );
  df <- df[ rows, ];

  if( nrow( df ) > 0 ) {
#    rn <- range(  df[[ xCountVar ]][ rows ], df[[ yCountVar ]][ rows ] );
    df$facet <- df[[ subjectVar ]];

    if( !is.null( colorVar ) ) {
      p <- ggplot( df, aes_string(
        x = xCountVar,
        y = yCountVar,
        color = if( !is.null( colorVar ) ) colorVar else NULL
      ) ) +
        geom_point( size = 1, alpha = 0.4 );
    } else {
      p <- ggplot( df, aes_string(
        x = xCountVar,
        y = yCountVar
      ) ) +
        geom_point( size = 1, alpha = 0.4, color = "grey" );
    }
    p <- p +
      geom_abline( slope = 1, color = "darkgray", linetype = 3 );
    if( showFit ) {
      p <- p +
        stat_smooth( method = lm, color = "black", linetype = 2, na.rm = TRUE );
    }
    if( countScale == "log10" ) {
      p <- p +
        scale_x_log10() + scale_y_log10();
    #  scale_x_log10( limits = rn ) + scale_y_log10( limits = rn );
    } else {
#      p <- p +
      #xlim( rn ) + ylim( rn );
    }
    p <- p +
      theme_bw() +
      xlab( paste0( countVar, " (", groupVar, " == ", groupLevels[[1]], ")" )) +
      ylab( paste0( countVar, " (", groupVar, " == ", groupLevels[[2]], ")" )) +
      facet_wrap( ~ facet );
    if( !is.null( colorVar ) ) {
      cg <- gsub( ".x$|.y$", "", colorVar );
      nm <- groups[[ cg ]]$name;
      if( length( grep( ".x$", colorVar ) ) > 0 ) nm <- paste0( nm, " (", groupVar, " == ", groupLevels[[1]], ")" );
      if( length( grep( ".y$", colorVar ) ) > 0 ) nm <- paste0( nm, " (", groupVar, " == ", groupLevels[[2]], ")" );
      if( cg %in% names( groups ) ) {
        colorValues <- groups[[ cg ]]$df$color;
        names( colorValues ) <- groups[[ cg ]]$df$level;
        p <- p + scale_color_manual( name = nm, values = colorValues );
      }
    }
    p
  } else {
    NULL
  }
}

## -----------------------------------------------------------------------
## SASCExp
## -----------------------------------------------------------------------

sascExp <- function(
  counts, sampleSheet,
  sampleVar = "sample", featureVar = "feature",
  mergeSamplesVars = NULL,
  ...
) {
  # ----- sanity checks -----
  stopifnot( is.data.frame( counts ) );
  stopifnot( featureVar %in% colnames( counts ) );

  if( is.null( sampleSheet ) ) {
    sampleSheet <- data.frame( sample = setdiff( colnames( counts ), featureVar ) );
    colnames( sampleSheet ) <- c( sampleVar );
  }
  stopifnot( sampleVar %in% colnames( sampleSheet ) );

  # ----- merging columns, if requested -----
  if( !is.null( mergeSamplesVars ) ) {
    ss <- sampleSheet %>% select_( .dots = c( mergeSamplesVars, sampleVar ) );
    ss[[ ".newSample" ]] <- factor( gsub( " ", "", apply( ss[ mergeSamplesVars ], 1, paste, collapse = "." ) ) );

    d <- counts %>%
      gather_( sampleVar, "rawCount", setdiff( colnames( counts ), featureVar ) ) %>%
      left_join( ss %>% select_( .dots = c( sampleVar, ".newSample" ) ), by = sampleVar );
    d <- d %>%
      group_by_( .dots = c( featureVar, ".newSample" ) ) %>%
      summarize( rawCount = sum( rawCount ) ) %>%
      ungroup();
    d <- d %>%
      spread_( ".newSample", "rawCount" );

    counts <- d;
    ss[[ sampleVar ]] <- ss[[ ".newSample" ]];
    sampleSheet <- ss %>% dplyr::select( -.newSample ) %>% distinct();
  }

  # -----
  groupVars <- c( setdiff( colnames( sampleSheet ), sampleVar ), sampleVar );
  rownames( sampleSheet ) <- sampleSheet[[ sampleVar ]];

  countsSampleNames <- setdiff( colnames( counts ), featureVar );
  sampleSheetSampleNames <- as.character( sampleSheet[[ sampleVar ]] );
  if( length( setdiff( sampleSheetSampleNames, countsSampleNames ) ) != 0 ) {
    stop( paste0(
      "Sample sheet specifies samples not present in the counts table: ",
      paste(
        setdiff( sampleSheetSampleNames, countsSampleNames ),
        collapse = ", "
      )
    ) );
  }

  # ----- handling summary data encoded in the count table -----
  rows <- grepl( "^__", counts[[ featureVar ]] );

  dcs <- data.frame(
    `total_count` = apply( counts[ , sampleSheetSampleNames ], 2, function( v ) sum( as.numeric( v ) ) )
  );
  dcs[[ sampleVar ]] <- rownames( dcs );

  m <- as.matrix( counts[ !rows, sampleSheetSampleNames ] );
  rownames( m ) <- counts[[ featureVar ]][ !rows ];
  names( dimnames( m ) ) <- c( featureVar, sampleVar );

  das <- data.frame(
    `on_feature` = apply( m[ , sampleSheetSampleNames ], 2, function( v ) sum( as.numeric( v ) ) )
  );
  das[[ sampleVar ]] <- rownames( das );

  m <- m[ rowSums( abs( m ) ) > 1, ];

  groupNames <- colnames( sampleSheet );

  if( sum( rows ) > 0 ) {
    s <- counts[ rows, c( featureVar, sampleSheetSampleNames ) ] %>%
      melt( id.vars = featureVar, variable.name = sampleVar, value.name = "count" );
    s[[ featureVar ]] <- gsub( "__", "", s[[ featureVar ]] );
    s <- s %>% dcast( paste0( sampleVar, " ~ ", featureVar ), value.var = "count" );
  } else {
    s <- NULL;
  }

  facFun <- function( df, n, lvs ) {
    if( !is.factor( df[[ n ]] ) ) {
      df[[ n ]] <- factor( df[[ n ]], levels = lvs );
    }
    df
  }

  sampleSheet <- facFun( sampleSheet, sampleVar, sampleSheetSampleNames );
  dcs <- facFun( dcs, sampleVar, levels( sampleSheet[[ sampleVar ]] ) );
  das <- facFun( das, sampleVar, levels( sampleSheet[[ sampleVar ]] ) );

  as <- sampleSheet %>% select_( .dots = sampleVar ) %>%
    left_join( dcs, by = sampleVar ) %>%
    left_join( das, by = sampleVar );
  if( !is.null( s ) ) {
    s <- facFun( s, sampleVar, levels( sampleSheet[[ sampleVar ]] ) );
    as <- as %>% left_join( s, by = sampleVar );
  }
  as <- facFun( as, sampleVar, levels( sampleSheet[[ sampleVar ]] ) ); # TODO: why is this needed?

  rawCountsExp( m, sampleSheet %>% droplevels(), groupNames = groupNames, alignmentSummary = as, ... );
}

## -----------------------------------------------------------------------
## RawCountsExp
## -----------------------------------------------------------------------

# TODO: sasc/rawExp should stop if duplicates found in feature names

rawCountsExp <- function( rawCounts, modelFrame, modelFormula = NULL,
  modelContrasts = NULL, minMeanLogCpm = NULL,
  groupNames = colnames( modelFrame ), alignmentSummary = NULL
) {
  stopifnot( is.data.frame( rawCounts ) || is.matrix( rawCounts ) );
  stopifnot( is.data.frame( modelFrame ) );
  stopifnot( colnames( rawCounts ) == rownames( modelFrame ) );

  lapply( colnames( modelFrame ), function( nm ) {
    v <- modelFrame[ nm ];
    if( is.factor( v ) ) {
      if( "" %in% levels( v ) ) stop( paste0( "Factor '", nm, " has a level which is empty string." ) );
    }
  } );

  badRows <- grep( "^__", rownames( rawCounts ), perl = TRUE );
  if( length( badRows ) > 0 ) {
    warning( paste0( "Removing rows with __: ", paste( rownames( rawCounts )[ badRows ], collapse = ", " ) ) );
    rawCounts <- rawCounts[ -badRows, ];
  }

  rawCounts <- as.matrix( rawCounts );

  if( !is.null( minMeanLogCpm ) ) {
    warning( "rawCountsExp: minMeanLogCpm is depreciated; use: x <- selectFeatures( x, minMeanLog2CpmFeatures( x, minMeanLog2Cpm = 2 ) );" );

    library( edgeR );
    edger <- calcNormFactors( DGEList( counts = rawCounts ), norm.method = "tmm", iteration = F );
    cnts <- cpm( edger, log = TRUE, prior.count = 2 );
    mCnts <- apply( cnts, 1, mean );
    rawCounts <- rawCounts[ mCnts >= minMeanLogCpm, ];
  }

  if( is.null( names( dimnames( rawCounts ) ) ) ) {
    names( dimnames( rawCounts ) ) <- c( "tag", "sample" );
  }

  rce <- list(
    # ----- central matrix -----
    rawCounts = rawCounts,

    # ----- cols annotations -----
    frame = modelFrame,
    design = NULL,
    alignmentSummary = alignmentSummary, # changes when rows are removed

    # ----- rows annotations -----
    tagsAnnotations = NULL,

    # ----- object annotations -----
    formula = modelFormula,
    contrasts = NULL,

    groups = list(),

    countVar = "rawCount",
    sampleVar = names( dimnames( rawCounts ) )[[2]],
    tagVar = names( dimnames( rawCounts ) )[[1]],
    subjectVar = NULL,
    groupVar = NULL
  )
  class( rce ) <- c( "RawCountsExp", "CountsExp", class( rce ) );

  for( gn in groupNames ) rce <- addGroup( rce, modelFrame, gn );
  if( !is.null( modelFormula ) ) {
    rce <- setModelFormula( rce, modelFormula, modelContrasts = modelContrasts );
  }

  rce
}

## -----------------------------------------------------------------------

getRawCountsExp.RawCountsExp <- function( x ) x;
getRawCounts.RawCountsExp <- function( x ) x$rawCounts;
getCounts.RawCountsExp <- function( x ) getRawCounts( x );
getGroups.RawCountsExp <- function( x ) x$groups;
getSampleSheet.RawCountsExp <- function( x ) x$frame;
getFormula.RawCountsExp <- function( x ) x$formula;
getDesign.RawCountsExp <- function( x ) x$design;
getContrasts.RawCountsExp <- function( x ) x$contrasts;
getTagsAnnotations.RawCountsExp <- function( x ) x$tagsAnnotations;
getSubjectVar.RawCountsExp <- function( x ) x$subjectVar;
getGroupVar.RawCountsExp <- function( x ) x$groupVar;
getSampleVar.RawCountsExp <- function( x ) x$sampleVar;
getTagVar.RawCountsExp <- function( x ) x$tagVar;
getCountVar.RawCountsExp <- function( x ) x$countVar;
getAlignmentSummary.RawCountsExp <- function( x ) x$alignmentSummary;

## -----------------------------------------------------------------------

setModelFormula <- function( rce, modelFormula, ..., modelContrasts = NULL ) {
  stopifnot( inherits( modelFormula, "formula" ) );

  rce$formula <- modelFormula;

  rce$design <- model.matrix( modelFormula, data = rce$frame );
  if( colnames( rce$design )[[1]] == "(Intercept)" ) colnames( rce$design )[[1]] <- "Intercept";
  colnames( rce$design ) <- make.names( colnames( rce$design ) );

  l <- as.list( substitute( list( ... ) ) )[ -1L ];
  if( !is.null( modelContrasts ) ) {
    rce$contrasts <- makeContrasts( contrasts = modelContrasts, levels = rce$design );
  } else if( length( l ) > 0 && !is.null( l[[1]] ) ) {
    rce$contrasts <- makeContrasts( ..., levels = rce$design );
  } else {
    l <- as.list( setdiff( colnames( rce$design )[-1], "Intercept" ) );
    l$levels <- rce$design;
    rce$contrasts <- do.call( makeContrasts, l );
  }
  rce
}

setSubjectGroupVars <- function( rce, subjectVar, groupVar ) {
  rce$subjectVar <- subjectVar;
  rce$groupVar <- groupVar;
  rce
}

## -----------------------------------------------------------------------

library( colorspace )
# For all sample sheet columns it creates uniform representation (assigned color, linetype, shape)
# Allows manual redefinition
# First table always to be merged by sampleId.
# Columns of the table (sampleId, linetype, color, shape)
addGroup <- function(
  aRawCountsExp, sampleSheet, id,
  name = NULL, levels. = NULL, labels = NULL, colors = NULL, linetypes = NULL,
  colorFun = function( n ) rainbow_hcl( n )
) {
  if( id %in% names( aRawCountsExp$groups ) ) {
    group <- aRawCountsExp$groups[[ id ]];
  } else {
    if( is.null( levels. ) ) {
      if( is.factor( sampleSheet[ ,id ] ) ) {
        levels. <- factor( levels( sampleSheet[ ,id ] ), levels( sampleSheet[ ,id ] ) );
      } else {
        levels. <- levels( as.factor( sampleSheet[ ,id ] ) );
      }
    }
    group <- list(
      id = id,
      name = id,
      df = data.frame(
        level = levels.,
        stringsAsFactors = FALSE
      )
    );
    rownames( group$df ) <- as.character( levels. );
    group$df$label <- as.character( levels. );
  }
  if( !is.null( name ) ) group$name <- name;

  if( !is.null( colors ) ) {
    if( is.null( levels. ) ) {
      group$df$color <- colors;
    } else {
      group$df$color[ levels. ] <- colors;
    }
  } else if( is.null( group$df$color ) ) {
    group$df$color <- colorFun( nrow( group$df ) );
  }

  if( !is.null( linetypes ) ) {
    if( is.null( levels. ) ) {
      group$df$linetype <- linetypes;
    } else {
      group$df$linetype[ levels. ] <- linetypes;
    }
  } else if( is.null( group$df$linetype ) ) {
    group$df$linetype <- rep_len( 1:6, nrow( group$df ) );
  }

  group$sf <- group$df[ as.character( sampleSheet[ ,id ] ), ];
  rownames( group$sf ) <- rownames( sampleSheet );

  aRawCountsExp$groups[[ id ]] <- group;
  aRawCountsExp
}

## -----------------------------------------------------------------------

selectFeatures.RawCountsExp <- function( x, features, toRemove = FALSE ) {
  rc <- x$rawCounts;
  keepRows <- rownames( rc ) %in% features;
  if( toRemove ) keepRows <- !keepRows;

  x$rawCounts <- rc[ keepRows,, drop = FALSE ];

  if( !is.null( x$tagsAnnotations ) ) {
    x$tagsAnnotations <- x$tagsAnnotations[ keepRows,, drop = FALSE ];
  }

  if( !is.null( x$alignmentSummary ) ) {
    stopifnot( colnames( rc ) == as.character( x$alignmentSummary[[ getSampleVar( x ) ]] ) );
    if( !( "on_removed_features" %in% colnames( x$alignmentSummary ) ) ) {
      x$alignmentSummary$on_removed_features <- 0;
    }

    orf <- apply( rc[ !keepRows,, drop = FALSE ], 2, function( v ) sum( as.numeric( v ) ) );
    x$alignmentSummary$on_removed_features <- x$alignmentSummary$on_removed_features + orf;
    x$alignmentSummary$on_feature <- x$alignmentSummary$on_feature - orf;
  }

  x
}

minMeanLog2CpmFeatures <- function( x, minMeanLog2Cpm = 2, prior.count = 2, norm.method = "tmm", lib.size = NULL, ... ) {
  rc <- getRawCounts( x );
  dge <- DGEList( counts = rc );
  edger <- calcNormFactors( dge, norm.method = norm.method, iteration = F, lib.size = lib.size, ... );
  cnts <- cpm( edger, log = TRUE, prior.count = prior.count );
  rownames( rc )[ apply( cnts, 1, mean ) >= minMeanLog2Cpm ];
}

minOccLog2CpmFeatures <- function( x, minOccFrac = 1/4, minLog2Cpm = 2, prior.count = 2, norm.method = "tmm", lib.size = NULL, ... ) {
  rc <- getRawCounts( x );
  dge <- DGEList( counts = rc );
  edger <- calcNormFactors( dge, norm.method = norm.method, iteration = F, lib.size = lib.size, ... );
  cnts <- cpm( edger, log = TRUE, prior.count = prior.count );
  rownames( rc )[ apply( cnts, 1, function( v ) sum( v >= minLog2Cpm ) ) >= minOccFrac * ncol( cnts ) ];
}

# OBSOLETE: will be removed 20161001
filterCountsExp <- function( rce, filterFun = function( rawCounts ) rawCounts ) {
  warning( "filterCountsExp is depreciated; use selectFeatures" );

  if( length( formalArgs( filterFun ) ) == 1 ) {
    warning( "filterCountsExp: filterFun with one arg has been obsoleted; add second arg for tagsAnnotations" );
    rce$rawCounts <- filterFun( rce$rawCounts );
  } else {
    rce$rawCounts <- filterFun( rce$rawCounts, rce$tagsAnnotations );
  }
  if( !is.null( rce$tagsAnnotations ) ) {
    rce$tagsAnnotations <- rce$tagsAnnotations[ rownames( rce$rawCounts ), ];
  }
  rce
}

## -----------------------------------------------------------------------

selectSamples <- function( x, samples, toRemove = FALSE ) {
  stop( "Not implemented yet." )
}

filterSamplesExp <- function( rce, filterFun = function( modelFrame ) modelFrame ) {
  rce$frame <- droplevels( filterFun( rce$frame ) );
  nm <- rownames( rce$frame );
  rce$rawCounts <- rce$rawCounts[ , nm ];
  rce$groups <- lapply( rce$groups, function( g ) {
    g$sf <- g$sf[ nm, ];
    g
  } );
  if( !is.null( rce$formula ) ) {
    #    rce <- setModelFormula( rce, modelFormula = rce$formula, modelContrasts = rce$contrasts );
    warning( "Temporary: Resetting model formula to empty; set it once again.")
  }
  rce
}

## -----------------------------------------------------------------------

countsData.RawCountsExp <- function( x, addSampleSheet = FALSE, addDesign = FALSE, addAlignmentSummary = FALSE, addTagsAnnotations = FALSE, ... ) {
  tagVar <- getTagVar( x );
  sampleVar <- getSampleVar( x );
  countVar <- getCountVar( x );

  sampleSheet <- getSampleSheet( x );
  stopifnot( is.factor( sampleSheet[[ sampleVar ]] ) ); # sanity

  d <- getCounts( x ) %>% as.data.frame %>% add_rownames( var = tagVar );
  d <- d %>% gather_( sampleVar, countVar, setdiff( colnames( d ), tagVar ) );
  d[[ tagVar ]] <- factor( d[[ tagVar ]] ); # TODO: order should come from getTags or getTagsAnn
  d[[ sampleVar ]] <- factor( d[[ sampleVar ]], levels = levels( sampleSheet[[ sampleVar ]] ) );

  if( addSampleSheet ) {
    d <- d %>% left_join( sampleSheet, by = c( sampleVar ) );
  }

  if( addDesign ) {
    dd <- getDesign( x ) %>% as.data.frame() %>% add_rownames( var = sampleVar ) %>% dplyr::select( -Intercept );
    cns <- setdiff( colnames( dd ), colnames( getSampleSheet( x ) ) );
    for( cn in cns ) dd[[ cn ]] <- factor( dd[[ cn ]] );
    d <- d %>% left_join( dd, by = c( sampleVar ) );
  }

  if( addAlignmentSummary ) {
    as <- getAlignmentSummary( x );
    if( !is.null( as ) ) {
      d <- d %>% left_join( as, by = c( sampleVar ) );
    }
  }

  if( addTagsAnnotations ) {
    ta <- getTagsAnnotations( x );
    if( !is.null( ta ) ) {
      ta[[ tagVar ]] <- factor( as.character( ta[[ tagVar ]] ), levels = levels( d[[ tagVar ]] ) );
      d <- d %>% left_join( ta, by = c( tagVar ) );
    }
  }

  d
}

countsSummary <- function( rawCountsExp ) {
  rawCounts <- getCounts( rawCountsExp );
  list(
    numOfSamples = ncol( rawCounts ),
    numOfTags = nrow( rawCounts )
  )
}

samplesSummary <- function( rawCountsExp ) {
  warning( "samplesSummary is obsoleted now; use tagCountsSummary instead" );

  rawCounts <- getCounts( rawCountsExp );
  do.call( rbind, apply( rawCounts, 2, function( v ){ data.frame(
    totalCount = sum( as.numeric( v ), na.rm = TRUE ),
    numNA = sum( is.na( v ), na.rm = TRUE ),
    numZero = sum( v == 0 ),
    minCount = min( v, na.rm = TRUE ),
    medianCount = median( v, na.rm = TRUE ),
    maxCount = max( v, na.rm = TRUE )
  ) } ) );
}

tagCountsSummary <- function( x ) {
  countsDensityData( x ) %>%
    group_by_( .dots = getSampleVar( x ) ) %>%
    summarize(
      minTagCount = min( rawCount ),
      medianTagCount = median( rawCount ),
      maxTagCount = max( rawCount ),
      tagCountsSum = sum( as.numeric( rawCount ) ),
      naTagCountsNum = sum( is.na( rawCount ) ),
      zeroTagCountsNum = sum( rawCount == 0 )
    );
}

summary.RawCountsExp <- function( x ) {
  l <- list();
  l$dim <- dim( getCounts( x ) );
  l$alignmentSummary <- getAlignmentSummary( x );
  l$tagCountsSummary <- tagCountsSummary( x );

  ta <- getTagsAnnotations( x );
  if( !is.null( ta ) && !is.null( ta$chromosome_name ) ) {
    l[[ paste0( getTagVar( x ), "2chr" ) ]] <- table( ta$chromosome_name, useNA = "always" );
  }

  l
}

## -----------------------------------------------------------------------

addLogicalAnnotation <- function( x, varName, tags ) {
  tagVar <- getTagVar( x );

  ta <- x$tagsAnnotations;
  if( is.null( ta ) ) {
    ta <- data.frame( .tags = getTags( x ) );
    colnames( ta ) <- c( tagVar );
  }
  ta[[ varName ]] <- FALSE;
  ta[[ varName ]][ as.character( ta[[ tagVar ]] ) %in% tags ] <- TRUE;
  x$tagsAnnotations <- ta;

  x
}

addRefflatAnnotations <- function( x, refflatFile ) {
  tagVar <- getTagVar( x );
  d <- read.table( refflatFile, header = TRUE, sep = "\t" ) %>% as.tbl();
  d <- d %>% dplyr::select( -bStart, -bEnd );
  d <- d %>% group_by( gene, chr, strand );
  if( "geneStart" %in% colnames( d ) ) {
    d <- d %>%
      summarize(
        meanExonsNum = mean( exonsNum ),
        meanSumExonsLengths = mean( sumExonLengths ),
        geneStart = min( geneStart ),
        geneEnd = max( geneEnd )
      );
  } else {
    d <- d %>%
      summarize(
        meanExonsNum = mean( exonsNum ),
        meanSumExonsLengths = mean( sumExonLengths )
      );
  }
  d <- d %>%
    ungroup() %>%
    mutate( gene = as.character( gene ) );

  tags <- as.character( getTags( x ) );
  dd <- data.frame( .id = tags, stringsAsFactors = FALSE );
  dd <- dd %>% left_join( d, by = c( ".id" = "gene" ) );

  dd <- dd[ match( tags, dd$.id ), ];
  colnames( dd )[ which( colnames( dd ) == ".id" ) ] <- tagVar;

  stopifnot( dd[[ tagVar ]] == getTags( x ) );

  x$tagsAnnotations <- dd;
  x
}

# ---

addEnsemblAnnotations <- function( x, refId, tagVar = "feature" ) {
  library( biomaRt );

  tagsNames <- as.character( getTags( x ) );

  if( refId == "human" ) {
    if( length( grep( "ENSG", tagsNames ) ) > 20 ) {
      tagAttr <- "ensembl_gene_id";
    } else {
      tagAttr <- "hgnc_symbol";
    }
    cleanChromosomeNames <- c( 1:22, "X", "Y", "MT" );
    attrs <- c( "ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position" );
    host <- "www.ensembl.org";
    dataset <- "hsapiens_gene_ensembl";
    biomart <- "ensembl";
    ret <- addEnsemblGeneAnnotations( x );
  } else if( refId == "mouse" ) {
    tst <- c( 'Arnt2', 'Ascl2', 'C2cd5', 'Cab39', 'Cap2', 'Ccp110', 'Cep152', 'Chmp5', 'Cldn4', 'Clic6',
      'Col6a5', 'Cul3', 'Dag1', 'Dcdc2b', 'Dolk', 'Eci1', 'Esd', 'Exosc9', 'Fam129b', 'Fancf', 'Fhl3', 'Fnbp4',
      'Fubp3', 'Gar1', 'Gata6', 'Gmppa', 'Grhl3', 'Gria3', 'Gtf2ird1', 'Gtsf1l', 'Hadha', 'Hspe1', 'Ide', 'Inpp5a',
      'Itgav', 'Kdsr', 'Kifc1', 'Klhdc3', 'Kmt2c', 'Lhpp', 'Limk2', 'Man1a', 'March2', 'Mbd3', 'Mmrn2', 'Msantd2',
      'Mthfr', 'Nphp1', 'Ocel1', 'Pign', 'Pitpna', 'Plcxd1', 'Pmm2', 'Ppp1ca', 'Ptcd2', 'Ptprs', 'Rap1gap',
      'Rgs16', 'Rmdn3', 'Rnf6', 'Rspry1', 'Sall3', 'Slc25a33', 'Sord', 'Sox11', 'Spata5l1', 'Spop', 'Spp1',
      'Sra1', 'Stt3a', 'Sumo1', 'Syt12', 'Taz', 'Tbcb', 'Tmed8', 'Tmem209', 'Trmt13', 'Tsnax', 'Tssc1', 'Ubtd2',
      'Wdfy3', 'Wdr31', 'Wdr75', 'Zcchc6', 'Zfp119b', 'Zfp637', 'Zyg11b'
    );
    ovl <- intersect( as.character( getTags( x ) ), tst );
    if( length( ovl ) > length( tst ) / 5 ) {
      tagAttr <- "mgi_symbol";
    } else {
      tagAttr <- "ensembl_gene_id";
    }
    stop( "Unknown refId. Supported: human" );

    biomart <- "ensembl";
    host <- "www.ensembl.org";
    dataset <- "mmusculus_gene_ensembl";
    ensembl <- useMart( biomart = biomart, host = host, dataset = dataset );

    attrs <- c( "mgi_symbol", "ensembl_gene_id", "ensembl_exon_id", "chromosome_name", "start_position", "end_position" );
    genesD <- getBM( attributes = union( attrs, tagAttr ), filters = tagAttr, values = tagsNames, mart = ensembl );

    attrs <- c( "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end" );
    exonsD <- getBM( attributes = attrs, filters = "ensembl_exon_id", values = unique( genesD$ensembl_exon_id ), mart = ensembl );

    genesD %>% dplyr::select( -ensembl_exon_id );
  } else {
    stop( "Unknown refId. Supported: human" );
  }
  ret
}

# for mouse: ddd <- addEnsemblGeneAnnotations( ddd, dataset = "mmusculus_gene_ensembl", renameToTagAttr = "external_gene_name", attrs = c( "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position" ), cleanChromosomeNames = c( 1:19, "X", "Y", "MT" ) );
# ddd <- addEnsemblGeneAnnotations( ddd, dataset = "mmusculus_gene_ensembl", tagAttr = "mgi_symbol", attrs = c( "ensembl_gene_id", "mgi_symbol", "chromosome_name", "start_position", "end_position" ), cleanChromosomeNames = c( 1:19, "X", "Y", "MT" ) );
addEnsemblGeneAnnotations <- function( rawCountsExp,
                                       tagsNames = as.character( getTags( rawCountsExp ) ),
                                       tagAttr = "ensembl_gene_id",
                                       renameToTagAttr = NULL,
                                       cleanChromosomeNames = c( 1:22, "X", "Y", "MT" ),
                                       attrs = c( "ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position" ),
                                       host = "www.ensembl.org", #"may2009.archive.ensembl.org/biomart/martservice/", # "www.ensembl.org"
                                       dataset = "hsapiens_gene_ensembl",
                                       biomart = "ensembl" # "ENSEMBL_MART_ENSEMBL"
) {
  library( biomaRt );

  ensembl <- useMart( biomart = biomart, host = host, dataset = dataset );

  results <- getBM( attributes = union( attrs, tagAttr ), filters = tagAttr, values = tagsNames, mart = ensembl );
  if( !is.null( cleanChromosomeNames ) ) {
    rows <- ( results$chromosome_name %in% cleanChromosomeNames );
    results$chromosome_name[ !rows ] = NA;
    results$chromosome_name <- factor( results$chromosome_name, levels = cleanChromosomeNames );
  }
  #d <- merge( data.frame( tagsName = tagsNames ), results, by.x = "tagsName", by.y = tagAttr, all.x = TRUE, all.y = FALSE );
  d <- results[ match( tagsNames, results[[ tagAttr ]] ), ];

  if( !is.null( renameToTagAttr ) ) {
    newTagsNames <- make.names( d[ ,renameToTagAttr ], unique = TRUE );
    pos <- is.na( newTagsNames );
    newTagsNames[ pos ] <- tagsNames[ pos ];
    rownames( rawCountsExp$rawCounts ) <- newTagsNames;
    tagsNames <- newTagsNames;
  }

  rownames( d ) <- tagsNames;
  rawCountsExp$tagsAnnotations <- d[ tagsNames, ];

  rawCountsExp
}

addGTFGeneAnnotations <- function( x, gtfFile, genome = NA, cleanChromosomeNames = paste0( "chr", c( 1:22, "X", "Y", "MT" ) ) ) {
  library( rtracklayer );
  d <- import( gtfFile, format = "gtf", genome = genome, feature.type = c( "exon" ) );
  d <- as.data.frame( d ) %>%
    filter( type == "exon" ) %>%
    dplyr::rename( chromosome_name = seqnames, geneId = gene_id );
  if( !is.null( cleanChromosomeNames ) ) d <- d %>%
    filter( chromosome_name %in% cleanChromosomeNames );
  d <- d %>%
    mutate( exonLength = end - start + 1 ) %>%
    group_by( geneId, transcript_id, chromosome_name, strand ) %>%
    summarise( start = min( start ), end = max( end ), exonsNum = n(), exonsLength = sum( exonLength ) );
  rownames( genes ) <- genes$geneId;
  class( genes )

  stop( "Not finished yet" );
}

addGeneAnnotations <- function( x, ann ) {
  stopifnot( inherits( ann, "data.frame" ) );

  ta <- getTagsAnnotations( x );
  if( is.null( ta ) ) {
    counts <- getCounts( x );
    ta <- data.frame( row.names = rownames( counts ) );
  }
  m <- merge( ta, ann, by = "row.names", all.x = TRUE, all.y = FALSE, sort = FALSE );
  rownames( m ) <- m$Row.names;
  x$tagsAnnotations <- m[ rownames( ta ), ];
  x
}

## -----------------------------------------------------------------------

addPredictedSex <- function( x,
                             yChrNames = c( "Y", "chrY" ),
                             excludeChrNames = c( "MT", "chrMT" ),
                             minMedian = 10,
                             colName = "predictedSex"
) {
  rawCounts <- getRawCounts( x );
  chr <- getTagsAnnotations( x )$chromosome_name;
  yRows <- !is.na( chr ) & chr %in% yChrNames;
  autosomalRows <- !is.na( chr ) & !( chr %in% c( yChrNames, excludeChrNames ) );
  autosomalN <- apply( rawCounts[ autosomalRows, ], 2, quantile, 0.05 );
  yN <- apply( rawCounts[ yRows, ], 2, median );
  predictedSex <- ifelse( autosomalN >= minMedian, ifelse(
    yN < minMedian, "female", "male"
  ), NA );
  x$frame[[ colName ]] <- predictedSex;
  x <- addGroup( x, x$frame, colName );
  x
}

## -----------------------------------------------------------------------

alignmentSummaryScatterPlot <- function( x,
  xVar = "total_count", yVar = "on_feature / total_count", colorVar = NULL, labelVar = getSampleVar( x )
) {
  sampleVar <- getSampleVar( x );
  sampleNames <- rownames( getSampleSheet( x ) );

  d <- getAlignmentSummary( x ) %>%
    left_join( getSampleSheet( x ), by = sampleVar );

  p <- ggplot( d, aes_string( x = xVar, y = yVar, color = colorVar, label = labelVar ) ) +
    theme_bw();
  if( is.null( labelVar ) ) {
    p <- p + geom_point();
  } else {
    p <- p + geom_text();
  }
  p
}

alignmentSummaryData <- function( x, groupVar = "readGroup", countVar = "readCount" ) {
  sampleVar <- getSampleVar( x );
  ss <- getSampleSheet( x );
  sampleNames <- ss[[ sampleVar ]];
  if( is.null( sampleNames ) ) sampleNames <- rownames( ss );

  d <- getAlignmentSummary( x ) %>% left_join( ss ) %>%
    melt( id.vars = unique( c( sampleVar, "total_count", colnames( ss ) ) ), variable.name = groupVar, value.name = countVar );
  d[[ sampleVar ]] <- factor( d[[ sampleVar ]], rev( sampleNames ) );

  d
}

alignmentSummaryPlot <- function( x, showFractions = FALSE,
                                  groupVar = "readGroup", countVar = "readCount",
                                  sampleVar = getSampleVar( x ), flip = TRUE
) {
  d <- alignmentSummaryData( x, groupVar = groupVar, countVar = countVar );

  p <- ggplot( d, aes_string( x = sampleVar, y = countVar, fill = groupVar ) );
  if( showFractions ) {
    p <- p + geom_bar( stat = "identity", position = "fill" ) +
      scale_y_continuous( labels = percent, name = "count fraction" );
  } else {
    p <- p + geom_bar( stat = "identity" );
  }
  p <- p + theme_bw();
  if( flip ) {
    p <- p + coord_flip();
  } else {
    p <- p + theme( axis.text.x = element_text( angle = -90, vjust = 0.5, hjust = 1 ) );
  }
  p
}

## -----------------------------------------------------------------------

complexityData <- function( x, rankVar = "Rank", fractionVar = "Fraction" ) {
  l <- list();
  l[[ fractionVar ]] <- ".cumFracCount";
  l[[ rankVar ]] <- ".rank";
  countsDensityData( getRawCountsExp( x ) ) %>%
    group_by_( .dots = c( getSampleVar( x ) ) ) %>%
    arrange( desc( rawCount ) ) %>%
    mutate( .cumFracCount = cumsum( rawCount ) / sum( rawCount ), .rank = 1:n() ) %>%
    rename_( .dots = l )
}

complexityPlot <- function( x, ..., maxRank = 1000, showLabels = FALSE, labelVar = getSampleVar( x ), colorVar = getSampleVar( x ), rankVar = "Rank", fractionVar = "Fraction" ) {
  d <- complexityData( x, rankVar = rankVar, fractionVar = fractionVar );
  d$asText <- FALSE;
  if( showLabels ) {
    d$asText <- unsplit(
      lapply(
        split( d, d$sample ),
        function( d ) {
          d$asText[ c( sample( 2:7, 1 ), sample( 8:25, 1 ), sample( 25:100, 1 ) ) ] <- TRUE;
          d$asText
        } ),
      d$sample
    );
  }
  d %>% filter( Rank <= maxRank );

  ggplot( d, aes_string( x = rankVar, y = fractionVar, color = colorVar, label = labelVar, group = labelVar ) ) +
    geom_line( alpha = 0.3 ) +
    geom_point( data = subset( d, asText == FALSE ), alpha = 0.3 ) +
    geom_text( data = subset( d, asText == TRUE ) ) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)),
      limits = c( 1, maxRank )
    ) +
    scale_y_continuous( labels = percent ) +
    ylab( "Cumulative fraction of total reads till Rank" ) +
    theme_bw() +
    annotation_logticks( sides = "tb" );
}

## -----------------------------------------------------------------------

frequentTagsData <- function( x, rankVar = "Rank", fractionVar = "Fraction", minFraction = 0.01, inAllSamples = TRUE ) {
  sampleVar <- getSampleVar( x );

  l <- list();
  l[[ fractionVar ]] <- ".cumFracCount";
  l[[ rankVar ]] <- ".rank";
  d <- countsDensityData( getRawCountsExp( x ) ) %>%
    group_by_( .dots = c( sampleVar ) ) %>%
    arrange( desc( rawCount ) ) %>%
    mutate( .cumFracCount = rawCount / sum( rawCount ), .rank = 1:n() ) %>%
    ungroup();
  if( inAllSamples ) {
    dd <- d %>%
      filter( .cumFracCount >= minFraction ) %>%
      dplyr::select_( .dots = c( getTagVar( x ) ) ) %>%
      distinct_( .dots = c( getTagVar( x ) ) );
    d <- d %>%
      semi_join( dd );
  } else {
    d <- d %>%
      filter( .cumFracCount >= minFraction );
  }
  d <- d %>% rename_( .dots = l );
  d[[ sampleVar ]] <- factor( d[[ sampleVar ]], as.character( getSampleSheet( x )[[ sampleVar ]] ) );
  d
}

frequentTagsPlot <- function( x, ..., minFraction = 0.01, rankVar = "Rank", fractionVar = "Fraction" , flip = FALSE ) {
  d <- frequentTagsData( x, rankVar = rankVar, fractionVar = fractionVar, minFraction = minFraction, ... );
  sampleVar <- getSampleVar( x );
  tagVar <- getTagVar( x );

  p <- ggplot( d, aes_string( x = sampleVar, y = fractionVar, fill = tagVar ) );
  p <- p + geom_bar( stat = "identity" );
  p <- p + geom_hline( yintercept = minFraction, color = "black", linetype = 2, alpha = 0.3 );
  p <- p + theme_bw();
  if( flip ) {
    p <- p + coord_flip();
  } else {
    p <- p + theme( axis.text.x = element_text( angle = -90, vjust = 0.5, hjust = 0 ) );
  }
  p
}

## -----------------------------------------------------------------------

countsBoxData <- function( x ) {
  tagCountsSummary( x ) %>% left_join( groupsSummary( x ) );
}

countsBoxPlot <- function( x,
  colorVar = .getGroupName( x, 1 ),
  coordVar = "tagCountsSum", colorGroup = NULL, labelVar = getSampleVar( x )
) {
  if( !is.null( colorGroup ) ) {
    warning( "boxplotCounts: colorGroup param is now called colorVar" );
    colorVar = colorGroup;
  }

  d <- countsBoxData( x );

  p <- ggplot( d, aes_string(
    y = coordVar,
    x = if( !is.null( colorVar ) ) colorVar else NULL,
    color = if( !is.null( colorVar ) ) colorVar else NULL,
    label = labelVar
  ) ) +
    geom_boxplot();
  if( !is.null( labelVar ) ) p <- p +
    geom_text( angle = 90, color = "black", size = 4 );

  if( 1 ) {
    p <- p + coord_trans( y = "log10" );
  } else {
    p <- p + scale_y_log10(
      breaks = trans_breaks( "log10", function(x) 10^x),
      labels = trans_format( "log10", math_format(10^.x) )
    );
  }
  p <- p +
    coord_flip() +
    theme_bw();

  groups <- getGroups( x );
  if( !is.null( colorVar ) ) {
    colorValues <- groups[[ colorVar ]]$df$color;
    names( colorValues ) <- groups[[ colorVar ]]$df$level;
    p <- p + scale_color_manual( name = groups[[ colorVar ]]$name, values = colorValues );
    p <- p + xlab( groups[[ colorVar ]]$name );
  }
  p
}

## -----------------------------------------------------------------------

countsDensityPlot.RawCountsExp <- function( x, ..., countScale = "log10" ) {
  countsDensityPlot.general( x, ..., countScale = countScale, zeroCount = 0.5 );
}

pairedCountsPlot.RawCountsExp <- function( x, ..., countScale = "log10" ) {
  pairedCountsPlot.general( x, ..., countScale = countScale, zeroCount = 0.5 );
}

## -----------------------------------------------------------------------

sexMixupPlot <- function( x,
  sexVar = "sex",
  yGenes = c( "NRY", "AZF1", "BPY2", "DAZ1", "DAZ2", "PRKY", "RBMY1A1", "SRY", "TSPY", "USP9Y", "UTY", "ZFY" ),
  countVar = getCountVar( x ),
  tagVar = getTagVar( x ),
  sampleVar = getSampleVar( x )
) {
  d <- countsData( x, addTagsAnnotations = TRUE, addDEPattern = FALSE, addDETab = FALSE );
  d <- d[ d[[ tagVar ]] %in% yGenes, ];
  d <- d %>% left_join( getSampleSheet( x ) );
  p <- ggplot( d, aes_string( y = sampleVar, x = countVar, color = tagVar ) ) +
    geom_point( alpha = 0.4 ) +
    theme_bw() +
    theme( axis.text.x = element_text( angle = -90, vjust = 0.5, hjust = 1 ) );
  p + facet_grid( paste0( ". ~ ", sexVar ) )
}

## -----------------------------------------------------------------------
## simpleNormalizedExp
## -----------------------------------------------------------------------

simpleNormalizedExp <- function( x, pseudoCount = 0.1, log = TRUE, ... ) {
  ret <- list(
    rawCountsExp = x,
    pseudoCount = pseudoCount,
    log = log
  );
  class( ret ) <- c( "SimpleNormalizedExp", "CountsExp", class( ret ) );
  ret
}

getRawCountsExp.SimpleNormalizedExp <- function( x ) x$rawCountsExp;
getNormalizedExp.SimpleNormalizedExp <- function( x ) x;

getCounts.SimpleNormalizedExp <- function( x, ... ) {
  rc <- getRawCounts( x ) + x$pseudoCount;
  s <- apply( rc, 2, sum, na.rm = TRUE );
  d <- t( t( rc ) / s );
  if( x$log ) d <- base::log10( d );
  names( dimnames( d ) ) <- names( dimnames( getRawCounts( x ) ) );
  d
}
getCountVar.SimpleNormalizedExp <- function( x ) "fraction";

## -----------------------------------------------------------------------
## edgeR: TMMCountsExp, edgeRFit
## -----------------------------------------------------------------------

tmmNormalizedExp <- function( rawCountsExp, ... ) {
  library( edgeR );

  rawCounts <- getRawCounts( rawCountsExp );
  ret <- list(
    rawCountsExp = rawCountsExp,
    edger = calcNormFactors( DGEList( counts = rawCounts, genes = rownames( rawCounts ) ), ... )
  );
  class( ret ) <- c( "TMMNormalizedExp", "CountsExp", class( ret ) );
  ret
}

getRawCountsExp.TMMNormalizedExp <- function( x ) x$rawCountsExp;
getNormalizedExp.TMMNormalizedExp <- function( x ) x;
getEdgeR.TMMNormalizedExp <- function( x ) x$edger;

getCounts.TMMNormalizedExp <- function( x, log = TRUE, prior.count = 2, ... ) {
  d <- cpm( x$edger, log = log, prior.count = prior.count, ... );
  names( dimnames( d ) ) <- names( dimnames( getRawCounts( x ) ) );
  d
}
getCountVar.TMMNormalizedExp <- function( x ) "log2CPM";

if( 0 ) { # obsoleted 20160614
  countsData.TMMNormalizedExp <- function( x, addTagsAnnotations = TRUE, rpkmVar = "log2RPKM", geneLenVar = "meanSumExonsLengths", log = TRUE, prior.count = 2, ... ) {
    tagVar <- getTagVar( x );
    sampleVar <- getSampleVar( x );

    d <- countsData( getRawCountsExp( x ), addTagsAnnotations = addTagsAnnotations, ... );

    f <- function( ddd, countVar ) {
      dd <- ddd %>% as.data.frame %>% add_rownames( var = tagVar );
      dd <- dd %>% gather_( sampleVar, countVar, setdiff( colnames( dd ), tagVar ) );
      dd[[ tagVar ]] <- factor( dd[[ tagVar ]], levels = levels( d[[ tagVar ]] ) );
      dd[[ sampleVar ]] <- factor( dd[[ sampleVar ]], levels = levels( d[[ sampleVar ]] ) );
      dd
    }

    ta <- getTagsAnnotations( x );

    cpmD <- cpm( x$edger, log = log, prior.count = prior.count );
    tags <- rownames( cpmD );
    ta <- ta[ match( ta[[ tagVar ]], tags ), ];
    stopifnot( ta[[ tagVar ]] == tags );
    d <- d %>% left_join( f( cpmD, getCountVar( x ) ), by = c( tagVar, sampleVar ) );

    if( !is.null( ta[[ geneLenVar ]] )) {
      rpkmD <- rpkm( x$edger, gene.length = ta[[ geneLenVar ]], log = log, prior.count = prior.count )
      stopifnot( rownames( rpkmD ) == tags );

      d <- d %>% left_join( f( rpkmD, rpkmVar ), by = c( tagVar, sampleVar ) );
    }

    d
  }
}

writeReport.TMMNormalizedExp <- function( x, names = c(), ... ) {
  tagVar <- getTagVar( x );
  sampleVar <- getSampleVar( x );

  f <- function( countVar, ... ) {
    d <- countsData( x );
    if( countVar %in% names( d ) ) {
      d <- d %>%
        select_( .dots = c( tagVar, sampleVar, countVar ) ) %>%
        spread_( sampleVar, countVar );

      writeReport( d, names = c( names, report = countVar ), ... );
    }
  }

  f( "log2CPM", ... );
  f( "log2RPKM", ... );
}

## -----------------------------------------------------------------------

edgeRFit.TMMNormalizedExp <- function( normalizedExp, df = 5, useRobust = FALSE, ... ) {
  edgeR <- getEdgeR( normalizedExp );
  design <- getDesign( normalizedExp );

  if( !useRobust ) {
    edgeR <- estimateGLMCommonDisp( edgeR, design );
    edgeR <- estimateGLMTrendedDisp( edgeR, design, df = df );
    edgeR <- estimateGLMTagwiseDisp( edgeR, design );
  } else {
    edgeR <- estimateGLMRobustDisp( edgeR, prior.df = df, ... );
  }
  fit <- glmFit( edgeR, design );

  ret <- list(
    normalizedExp = normalizedExp,
    fit = fit,
    edgeR = edgeR,
    useRobust = useRobust
  );
  class( ret ) <- c( "EdgeRFit", "CountsExp", class( ret ) );
  ret
}

getRawCountsExp.EdgeRFit <- function( x ) getRawCountsExp( x$normalizedExp );
getNormalizedExp.EdgeRFit <- function( x ) getNormalizedExp( x$normalizedExp );
getCounts.EdgeRFit <- function( x ) getCounts( x$normalizedExp );
getCountVar.EdgeRFit <- function( x ) getCountVar( x$normalizedExp );

## -----------------------------------------------------------------------

deTable.EdgeRFit <- function(
  x,
  contrastNames = colnames( getContrasts( x ) ),
  threshold = 0.05,
  adjust.method = "BH"
) {
  contrasts <- getContrasts( x );
  d <- do.call( rbind, lapply( contrastNames, function( contrastName ) {
    d <- glmLRT( x$fit, contrast = contrasts[ , contrastName ] );
    stopifnot( rownames( d$table ) == as.character( getTags( x ) ) );
    d <- with( d$table, data.frame( avgLog2CPM = logCPM, avgLog2FC = logFC, pValue = PValue ) );
    d[[ getTagVar( x ) ]] <- getTags( x );
    d$adj.P.Val <- p.adjust( d$pValue, method = adjust.method );
    d$DE <- factor( ( d$adj.P.Val <= threshold ) * sign( d$avgLog2FC ), levels = c( -1, 0, 1 ) );
    d$Rank <- NA_integer_;
    d$Rank[ order( d$adj.P.Val, d$pValue ) ] <- 1:nrow( d );
    d$Contrast <- contrastName;
    d
  } ) );

  ret <- deTable.general( x, deTab = d, contrastNames = contrastNames, threshold = threshold );
  class( ret ) <- c( "EdgeRDETable", class( ret ) );
  ret
}

plotDiag.EdgeRDETable <- function( x, ... ) {
  plotBCV( x$fit$edgeR );
}

getRawCountsExp.EdgeRDETable <- function( x ) getRawCountsExp( x$fit );
getNormalizedExp.EdgeRDETable <- function( x ) getNormalizedExp( x$fit );
getCounts.EdgeRDETable <- function( x ) getCounts( x$fit );
getCountVar.EdgeRDETable <- function( x ) getCountVar( x$fit );
getAvgCountVar.EdgeRDETable <- function( x ) "avgLog2CPM";
getDiffCountVar.EdgeRDETable <- function( x ) "avgLog2FC";

## -----------------------------------------------------------------------

edgeRAnalysis <- function( x, threshold = 0.05, ... ) {
  nce <- tmmNormalizedExp( x );
  fit <- edgeRFit( nce, ... );
  deTable( fit, threshold = threshold );
}

## -----------------------------------------------------------------------
## Voom
## -----------------------------------------------------------------------

voomNormalizedExp <- function( x, useTmmNormFactors = TRUE, ... ) {
  if( useTmmNormFactors ) {
    rc <- getRawCounts( x );
    edger <- calcNormFactors( DGEList( counts = rc ), norm.method = "tmm", iteration = F );
    v <- voom(
      counts = rc,
      design = getDesign( x ),
      plot = FALSE,
      lib.size = colSums( rc ) * edger$samples$norm.factors,
      ...
    );
  } else {
    v <- voom(
      counts = getRawCounts( x ),
      design = getDesign( x ),
      plot = FALSE,
      ...
    )
  }
  ret <- list(
    rawCountsExp = x,
    voom = v
  );
  class( ret ) <- c( "VoomNormalizedExp", "CountsExp", class( ret ) );
  ret
}

getRawCountsExp.VoomNormalizedExp <- function( x ) getRawCountsExp( x$rawCountsExp );
getNormalizedExp.VoomNormalizedExp <- function( x ) x;
getVoom.VoomNormalizedExp <- function( x ) x$voom;
getCounts.VoomNormalizedExp <- function( x ) x$voom$E;
getCountVar.VoomNormalizedExp <- function( x ) "log2CPM";
getWeights.VoomNormalizedExp <- function( x ) {
  w <- x$voom$weights;
  colnames( w ) <- colnames( getCounts( x ) );
  rownames( w ) <- rownames( getCounts( x ) );
  names( dimnames( w ) ) <- names( dimnames( getCounts( x ) ) );
  w
}

## -----------------------------------------------------------------------

eBayesFit.VoomNormalizedExp <- function( normalizedExp ) {
  ret <- list(
    normalizedExp = normalizedExp,
    fit = eBayes(
      contrasts.fit(
        lmFit( getVoom( normalizedExp ), getDesign( normalizedExp ) ),
        getContrasts( normalizedExp )
      )
    )
  );
  class( ret ) <- c( "EBayesFit", "CountsExp", class( ret ) );
  ret
}

getRawCountsExp.EBayesFit <- function( x ) getRawCountsExp( x$normalizedExp );
getNormalizedExp.EBayesFit <- function( x ) getNormalizedExp( x$normalizedExp );
getCounts.EBayesFit <- function( x ) getCounts( x$normalizedExp );
getCountVar.EBayesFit <- function( x ) getCountVar( x$normalizedExp );

## -----------------------------------------------------------------------

deTable.EBayesFit <- function( fit, contrastNames = colnames( getContrasts( fit ) ), threshold = 0.05 ) {
#  contrasts <- getContrasts( fit );
  d <- do.call( rbind, lapply( contrastNames, function( contrastName ) {
    d <- topTable( fit$fit, coef = contrastName, number = +Inf, sort.by = "none", confint = TRUE );
    if( "ID" %in% colnames( d ) ) rownames( d ) <- d$ID;
    colnames( d )[ which( colnames( d ) == "logFC" ) ] <- "avgLog2FC";
    colnames( d )[ which( colnames( d ) == "AveExpr" ) ] <- "avgLog2CPM";
    stopifnot( rownames( d ) == as.character( getTags( fit$normalizedExp ) ) );
    d[[ getTagVar( fit ) ]] <- getTags( fit$normalizedExp );
    d$Rank <- NA_integer_;
    d$Rank[ order( d$adj.P.Val, d$P.Value )] <- 1:nrow( d );
    d$DE <- factor( ( d$adj.P.Val <= threshold ) * sign( d$avgLog2FC ), levels = c( -1, 0, 1 ) );
    d$Contrast <- contrastName;
    d
  } ) );

  ret <- deTable.general( fit, deTab = d, contrastNames = contrastNames, threshold = threshold );
  class( ret ) <- c( "EBayesDETable", class( ret ) );
  ret
}

getRawCountsExp.EBayesDETable <- function( x ) getRawCountsExp( x$fit );
getNormalizedExp.EBayesDETable <- function( x ) getNormalizedExp( x$fit );
getCounts.EBayesDETable <- function( x ) getCounts( x$fit );
getCountVar.EBayesDETable <- function( x ) getCountVar( x$fit );
getAvgCountVar.EBayesDETable <- function( x ) "avgLog2CPM";
getDiffCountVar.EBayesDETable <- function( x ) "avgLog2FC";

## -----------------------------------------------------------------------

voomAnalysis <- function( rce ) {
  nce <- voomNormalizedExp( rce );
  fit <- eBayesFit( nce );
  deTable( fit );
}

## -----------------------------------------------------------------------

gsaAll.VoomNormalizedExp <- function( x, contrastNames = colnames( getContrasts( x ) ) ) {
  library( globaltest );

  counts <- getCounts( x );
  contrasts <- getContrasts( x );
  design <- getDesign( x );
  nullDesign <- design[ , !( colnames( design ) %in% contrastNames ), drop = FALSE ];

  l <- lapply( contrastNames, function( contrastName ) {
    gt( design[ , contrastName ], t( counts ), null = nullDesign, standardize = TRUE );
  } );
  names( l ) <- contrastNames;
  l
}

gsaKEGG.VoomNormalizedExp <- function( x, contrastNames = colnames( getContrasts( x ) ),
                                       annotation = "org.Hs.eg.db",
                                       id_mapping = ( function(){ library( org.Hs.eg.db ); as.list( org.Hs.egSYMBOL2EG ) } )()
) {
  library( KEGG.db );

  counts <- getCounts( x );
  contrasts <- getContrasts( x );
  design <- getDesign( x );
  nullDesign <- design[ , !( colnames( design ) %in% contrastNames ), drop = FALSE ];

  #symbol_mapping <- as.list(setNames(getGeneNames(rownames(countdata.filt), id.type="ensembl"),
  #                                   rownames(countdata.filt)))

  l <- lapply( contrastNames, function( contrastName ) {
    gtKEGG( design[ , contrastName ], t( counts ), null = nullDesign, annotation = annotation, probe2entrez = id_mapping );
  } );
  names( l ) <- contrastNames;

  setsTab <- rbind_all( lapply( contrastNames, function( contrastName ) {
    d <- subset( result( l[[ contrastName ]] ), holm <= 0.05 );
    if( nrow( d ) > 0 ) {
      d$Contrast <- contrastName;
      d$Index <- 1:nrow( d );
      d$KEGG <- rownames( d );
    } else {
      d <- data.frame();
    }
    d
  } ) );
  if( nrow( setsTab ) == 0 ) setsTab <- data.frame( Contrast = factor( c(), levels = contrastNames ) );

  genesTab <- rbind_all( lapply( contrastNames, function( contrastName ) {
    d <- subset( setsTab, Contrast == contrastName );
    if( nrow( d ) > 0 ) {
      ids <- d$KEGG;
      ll <- l[[ contrastName ]][ ids, ];
      d <- rbind_all( lapply( 1:length( ll ), function( n ) {
        d <- result( covariates( ll[[n]], plot = FALSE ) );
        data.frame( d, tag = rownames( d ), Contrast = contrastName, geneset = names( ll )[[n]], gtSetRank = n, gtTagRank = 1:nrow( d ) )
      } ) );
    } else {
      d <- data.frame();
    }
    d
  } ) );
  if( nrow( genesTab ) == 0 ) genesTab <- data.frame( tag = c(), Contrast = factor( c(), levels = contrastNames ) );

  ret <- list(
    gts = l,
    setsTab = setsTab,
    genesTab = genesTab
  );
  ret;
}

gsaKEGG.EBayesDETable <- function( x, ... ) {
  ret <- gsaKEGG.VoomNormalizedExp( x, ... );
  contrastNames = colnames( getContrasts( x ) );

  ret$genesTab <- rbind_all( lapply( split( ret$genesTab, ret$genesTab$Contrast ), function( gt ) {
    if( nrow( gt ) > 0 ) {
      d <- merge( gt, deData( x ), by = c( "tag", "Contrast" ) );
      d[ order( d$gtSetRank, d$gtTagRank ), ];
    } else {
      data.frame();
    }
  } ) );
  if( nrow( ret$genesTab ) == 0 ) ret$genesTab <- data.frame( tag = c(), Contrast = factor( c(), levels = contrastNames ) );

  class( ret ) <- c( "GsaKEGG", class( ret ) );
  ret
}

writeReport.GsaKEGG <- function( x, names = c(), ... ) {
  writeReport( x$setsTab, names = c( names, gsa = "KEGG", report = "setsTab" ), ... );
  writeReport( x$genesTab, names = c( names, gsa = "KEGG", report = "genesTab" ), ... );
}

plot.GsaKEGG <- function( x, contrastNames = unique( x$setsTab$Contrast ), labelVar = "Rank", colorVar = NULL, ... ) {
  setsTab <- x$setsTab;
  ord <- paste0( setsTab$Contrast, ":", setsTab$KEGG );

  deTab <- subset( x$genesTab, Contrast %in% contrastNames );

  if( nrow( deTab ) > 0 ) {
    deTab$ord <- factor( paste0( deTab$Contrast, ":", deTab$geneset ), levels = ord );
    p <- ggplot( deTab, aes( x = log2AvgExpr, y = log2FC, ymin = CI.L, ymax = CI.R, color = log10( p.value ) ) ) +
      geom_hline( yintercept = 1, color = "blue", linetype = 2, alpha = 0.5 ) +
      geom_hline( yintercept = -1, color = "blue", linetype = 2, alpha = 0.5 ) +
      geom_hline( yintercept = 0, color = "darkgray", linetype = 2, alpha = 0.5 ) +
      geom_errorbar( data = subset( deTab, p.value > 0.05 ), alpha = 0.5, color = "gray" ) +
      geom_errorbar( data = subset( deTab, p.value <= 0.05 ), alpha = 0.6 ) +
      geom_point( data = subset( deTab, p.value <= 0.05 ) ) +
      facet_wrap( ~ ord ) +
      theme_bw();
    p
  } else {
    NULL
  }
}

## -----------------------------------------------------------------------
## DESeq2CountsExp
## -----------------------------------------------------------------------

#  dne <- deseqNormalizedExp( rce, method = "blind", sharingMode = "fit-only" );
deseq2NormalizedExp <- function( rawCountsExp, ... ) {
  if( 0 ) {
    deseq <- DESeqDataSetFromMatrix(
      countData = getRawCounts( rawCountsExp ),
      colData = getSampleSheet( rawCountsExp ),
      design = getFormula( rawCountsExp )
    );
  } else {
#    .diffexpr.DESeqDataSet <- function (se, design, modelMatrix = model.matrix(design, data = as.data.frame(colData(se))) ) {

    se <- SummarizedExperiment( SimpleList( counts = getRawCounts( rawCountsExp ) ) );
    colData( se ) <- DataFrame( getSampleSheet( rawCountsExp ) );
    design <- getFormula( rawCountsExp );

    mcolsCols <- DataFrame(type = rep("input", ncol(colData(se))), description = rep("", ncol(colData(se))))

    mcols(colData(se)) <- if (is.null(mcols(colData(se)))) {
      mcolsCols
    } else if (all(names(mcols(colData(se))) == c("type", "description"))) {
      mcolsCols
    } else {
      cbind(mcols(colData(se)), mcolsCols)
    }

    dds <- new("DESeqDataSet", se, design = design)

    mcolsRows <- DataFrame(type = rep("input", ncol(mcols(dds))),
                           description = rep("", ncol(mcols(dds))))
    mcols(mcols(dds)) <- if (is.null(mcols(mcols(dds)))) {
      mcolsRows
    } else if (all(names(mcols(mcols(dds))) == c("type", "description"))) {
      mcolsRows
    } else {
      cbind(mcols(mcols(dds)), mcolsRows)
    }

    deseq <- dds;
  }

  deseq <- estimateSizeFactors( deseq );
  deseq <- estimateDispersions( deseq, ... );

  ret <- list(
    rawCountsExp = rawCountsExp,
    deseq = deseq
  );
  class( ret ) <- c( "DESeq2NormalizedExp", "CountsExp", class( ret ) );
  ret
}

getRawCountsExp.DESeq2NormalizedExp <- function( x ) getRawCountsExp( x$rawCountsExp );
getNormalizedExp.DESeq2NormalizedExp <- function( x ) x;
getDESeq2.DESeq2NormalizedExp <- function( x ) x$deseq;
getCounts.DESeq2NormalizedExp <- function( x ) {
  d <- log2( counts( getDESeq2( x ), normalized = TRUE ) );
  names( dimnames( d ) ) <- names( dimnames( getRawCounts( x ) ) );
  d
}
getCountVar.DESeq2NormalizedExp <- function( x ) "log2CPM";

# This is wrong, another function needed (plotFit)
#plot.DESeq2NormalizedExp <- function( x, ... ) {
#  plotDispEsts( getDESeq2( x ) ); grid()
#}

## -----------------------------------------------------------------------

nBinomWaldFit.DESeq2NormalizedExp <- function( x, ... ) {
  ret <- list(
    normalizedExp = x,
    fit = nbinomWaldTest( getDESeq2( x ), betaPrior = FALSE, ... )
  );
  class( ret ) <- c( "NBinomWaldFit", "CountsExp", class( ret ) );
  ret
}

getRawCountsExp.NBinomWaldFit <- function( x ) getRawCountsExp( x$normalizedExp );
getNormalizedExp.NBinomWaldFit <- function( x ) getNormalizedExp( x$normalizedExp );
getCounts.NBinomWaldFit <- function( x ) getCounts( x$normalizedExp );
getCountVar.NBinomWaldFit <- function( x ) getCountVar( x$normalizedExp );

## -----------------------------------------------------------------------

deTable.NBinomWaldFit <- function( fit, contrastNames = colnames( getContrasts( fit ) ), threshold = 0.05 ) {
  contrasts <- getContrasts( fit );
  d <- do.call( rbind, lapply( contrastNames, function( contrastName ) {
    d <- results( fit$fit, contrasts[ , contrastName ] );
    stopifnot( rownames( d ) == as.character( getTags( fit$normalizedExp ) ) );
    d <- with( d, data.frame(
      log2AvgExpr = log2( baseMean ), log2FC = log2FoldChange, pValue = pvalue, adj.P.Val = padj
    ) );
    d[[ getTagVar( fit ) ]] <- getTags( fit$normalizedExp );
    rownames( d ) <- d$tag;
    d$Rank <- NA_integer_;
    d$Rank[ order( d$adj.P.Val, d$pValue ) ] <- 1:nrow( d );
    d$DE <- factor( ( d$adj.P.Val <= threshold ) * sign( d$log2FC ), levels = c( -1, 0, 1 ) );
    num <- sum( is.na( d$DE ) );
    if( num > 0 ) {
      warning( paste0( "deTable.NBinomWaldFit: DESeq2 led to ", num, " genes with adjusted p-value(s) == NA; treating them as non-significant DE" ) );
      d$DE[ is.na( d$DE ) ] <- 0;
    }
    d$Contrast <- contrastName;
    d
  } ) );

  ret <- deTable.general( fit, deTab = d, contrastNames = contrastNames, threshold = threshold );
  class( ret ) <- c( "NBinomWaldDETable", class( ret ) );
  ret
}

getRawCountsExp.NBinomWaldDETable <- function( x ) getRawCountsExp( x$fit );
getNormalizedExp.NBinomWaldDETable <- function( x ) getNormalizedExp( x$fit );
getCounts.NBinomWaldDETable <- function( x ) getCounts( x$fit );
getCountVar.NBinomWaldDETable <- function( x ) getCountVar( x$fit );

## -----------------------------------------------------------------------

deseq2Analysis <- function( rce ) {
  nce <- deseq2NormalizedExp( rce );
  fit <- nBinomWaldFit( nce );
  deTable( fit );
}

## -----------------------------------------------------------------------
## DETable
## -----------------------------------------------------------------------

deTable.general <- function( x, deTab, contrastNames, threshold ) {
  tagVar <- getTagVar( x );

  #  colNames <- c( tagVar, "Contrast", "log2AvgExpr", "log2FC", "pValue", "adj.P.Val", "DE", "Rank" );

  #  deTab$DE <- factor( ( deTab$adj.P.Val <= threshold ) * sign( deTab$log2FC ), levels = c( -1, 0, 1 ) );
  #  d$Rank <- NA_integer_;
  #  d$Rank[ order( d$Contrast, d$adj.P.Val, d$pValue )] <- 1:nrow( d );
  #  tagAnn <- getTagsAnnotations( x$normalizedExp );
  #  if( !is.null( tagAnn ) ) d <- data.frame( d, tagAnn );

  deTab$Contrast <- factor( deTab$Contrast, levels = contrastNames );

  ret <- list(
    fit = x,
    deTab = deTab,
    threshold = threshold
  );
  class( ret ) <- c( "DETable", "CountsExp", class( ret ) );
  ret
}

addDEPattern <- function( x,
  contrastNames = colnames( getContrasts( x ) ),
  dePatternVar = paste0( contrastNames, collapse = "." )
) {
  tagVar <- getTagVar( x );

  d <- x$deTab;
  levels( d$DE ) <- list( "0" = 0, "u" = +1, "d" = -1 );
  d <- d %>% dcast( paste0( tagVar, " ~ Contrast" ), value.var = "DE" );
  d %>% select_( .dots = c( tagVar, contrastNames ) );

  refLevel <- paste0( rep( "0", length( contrastNames ) ), collapse = "." );
  d[[ dePatternVar ]] <- apply( d[ ,contrastNames, drop = FALSE ], 1, function( v ){ paste0( v, collapse = "." ); } );
  levs <- unique( c( refLevel, d[[ dePatternVar ]] ) );
  d[[ dePatternVar ]] <- relevel( factor( d[[ dePatternVar ]], levels = levs ), ref = refLevel );

  x$deTab <- left_join( x$deTab, d[ c( tagVar, dePatternVar ) ], by = tagVar );
  x
}

# ------------------------------------------------------------------------

# deData returns one row per Tag*Contrast
# should support maxRank = 50, topTagsNum = maxRank, contrastNames = colnames( getContrasts( x ) ), DEs = c( 1, -1 )
deData.DETable <- function( x, addTagsAnnotations = FALSE, geneLenVar = "meanSumExonsLengths", selectTags = NULL, ... ) {
  tagVar <- getTagVar( x );
  avgCountVar <- getAvgCountVar( x );

  d <- x$deTab;
  tags <- getTags( x );
  levels( d[[ tagVar ]] ) <- tags;

  stopifnot( is.null( selectTags ) );

  if( avgCountVar == "avgLog2CPM" && !is.na( geneLenVar ) ) {
    addTagsAnnotations <- TRUE;
  }

  if( addTagsAnnotations ) {
    ta <- getTagsAnnotations( x );
    if( !is.null( ta ) ) {
      #ta[[ tagVar ]] <- factor( as.character( ta[[ tagVar ]] ), levels = levels( d[[ tagVar ]] ) );
      d <- d %>% left_join( ta, by = c( tagVar ) );
    }
  }

  if( avgCountVar == "avgLog2CPM" && !is.na( geneLenVar ) && ( geneLenVar %in% colnames( d ) ) ) {
    d[[ "avgLog2RPKM" ]] <- d[[ avgCountVar ]] - log2( d[[ geneLenVar ]] / 1000 );
  }

  d
}

# countsData returns one row per Tag*Sample
countsData.DETable <- function( x,
  addDePattern = TRUE,
  dePatternVar = paste0( colnames( getContrasts( x ) ), collapse = "." ),
  addDeTab = TRUE,
  ...
) {
  tagVar <- getTagVar( x );
  sampleVar <- getSampleVar( x );
  countVar <- getCountVar( x );
  avgCountVar <- getAvgCountVar( x );

  d <- countsData( getNormalizedExp( x ), ... );

  if( addDeTab || addDePattern ) {
    ded <- deData( x );
    vars <- c( tagVar, avgCountVar );
    if( avgCountVar == "avgLog2CPM" && "avgLog2RPKM" %in% colnames( ded ) ) { # ugly, will be corrected
      vars <- c( vars, "avgLog2RPKM" );
    }
    if( dePatternVar %in% names( ded ) ) {
      vars <- c( vars, dePatternVar );
    }
    d <- d %>% left_join(
      ded %>% select_( .dots = vars ) %>% distinct()
    );
  }

  d
}

# ------------------------------------------------------------------------

# will be obsoleted soon -> deData
topTags.DETable <- function( x, maxRank = 50,
  topTagsNum = maxRank, contrastNames = colnames( getContrasts( x ) ), DEs = c( 1, -1 )
) {
  d <- deData( x );
  d <- subset( d, ( Contrast %in% contrastNames ) & ( DE %in% DEs ) & ( Rank <= topTagsNum ) );
  cols <- which( colnames( d ) %in% c( "t", "B", "ID" ) );
  if( length( cols ) >= 1 ) d <- d[ -cols ];
  d <- d[ order( d$Contrast, d$Rank ),, drop = FALSE ];
  rownames( d ) <- NULL;

  cols <- c( "Contrast", "Rank", getTagVar( x ), getAvgCountVar( x ), getDiffCountVar( x ), "DE", "CI.L", "CI.R" );
  cols <- intersect( cols, colnames( d ) );
  d <- data.frame( d[ cols ], d[ !( colnames( d ) %in% cols ) ] );
  d$Contrast <- factor( d$Contrast, levels = contrastNames );

  d
}

summary.DETable <- function( x, contrastNames = colnames( getContrasts( x ) ), ... ) {
  list(
    deContrasts = getContrasts( x ),
    deCounts = with( subset( deData( x ), Contrast %in% contrastNames ), ftable( Contrast, DE ) ),
    deGroupCounts = diffTags( x, ... )
  )
}

writeReport.DETable <- function( x, names = c(), ... ) {
  writeReport( deData( x ), names = c( names, report = "deTab" ), ... );
  writeReport( getNormalizedExp( x ), names = c( names, report = "normalizedExp" ), ... );
}

## -----------------------------------------------------------------------

diffTags <- function( x, contrastNames = colnames( getContrasts( x ) ), ... ) {
  tagVar <- getTagVar( x );
  d <- deData( x ) %>%
    filter( Contrast %in% contrastNames ) %>%
    dplyr::select_( .dots = c( tagVar, "Contrast", "DE" ) ) %>%
    spread( "Contrast", "DE" ) %>% # dcast( paste0( tagVar, " ~ Contrast" ), fill = 0, value.var = "DE" ) %>%
    dplyr::group_by_( .dots = contrastNames ) %>%
    summarize( count = n() ) %>% ungroup() %>%
    arrange( -count );
  d$Rank <- 1:nrow( d );
  d
}

contrastsDEData <- function( x, contrastNames = colnames( getContrasts( x ) ), removeNonDE = length( contrastNames ) > 1, dePatternVar = "patternVar", ... ) {
  tagVar <- getTagVar( x );
  d <- deData( x ) %>%
    filter( Contrast %in% contrastNames ) %>%
    dplyr::select_( .dots = c( tagVar, "Contrast", "DE" ) );
  if( removeNonDE ) {
    d <- d %>% dplyr::group_by_( .dots = tagVar ) %>%
    mutate( nonzeroDE = sum( DE != "0" ) ) %>%
    ungroup() %>%
    filter( nonzeroDE > 0 )
  }
  d <- d %>%
    dcast( paste0( tagVar, " ~ Contrast" ), fill = 0, value.var = "DE" );
  if( !is.null( dePatternVar ) ) {
    dd <- d[ contrastNames ];
    for( cn in contrastNames ) {
      dd[[ cn ]] <- factor( dd[[ cn ]] );
      levels( dd[[ cn ]] ) <- list( "d" = "-1", "0" = "0", "u" = 1 );
    }
    d[[ dePatternVar ]] <- apply( dd, 1, paste, collapse = ":" );
  }
  d
}

contrastsDECountData <- function( x, contrastNames = colnames( getContrasts( x ) ), removeNonDE = length( contrastNames ) > 1, dePatternVar = "patternVar", ... ) {
  d <- contrastsDEData( x, contrastNames = contrastNames, removeNonDE = removeNonDE, dePatternVar = dePatternVar, ... );

  d <- d %>%
    dplyr::group_by_( .dots = contrastNames ) %>%
    summarize( count = n() ) %>%
    ungroup() %>%
    arrange( desc( count ) );

  d$group <- 1:nrow( d );

  d <- d %>% melt( measure.vars = contrastNames, variable.name = "Contrast", value.name = "DE" );
  d$DE <- factor( d$DE, levels = c( "-1", "0", "1" ) );
  d
}

contrastsDECountPlot <- function( x, contrastNames = colnames( getContrasts( x ) ), ... ) {
  d <- contrastsDECountData( x, contrastNames = contrastNames, ... );
  ggplot( d, aes( x = Contrast, y = count, fill = DE, group = group ) ) +
    geom_bar( stat = "identity", color = "gray", fill = NA ) +
    geom_bar( stat = "identity", color = NA, alpha = 0.5 ) +
    theme_bw() +
    scale_fill_manual( values = c( "0" = NA, "1" = "red", "-1" = "blue" ) );
}

# ------------------------------------------------------------------------

# topTags.DETable->deData: consider which arguments are needed maxRank = 50, topTagsNum = maxRank, contrastNames = colnames( getContrasts( x ) ), DEs = c( 1, -1 )
chrDECountPlot <- function( x, topTagsNum = +Inf, ... ) {
  deTab <- topTags.DETable( x, topTagsNum = topTagsNum, DEs = c( -1, 0, 1, NA ), ... );
  deTab$chromosome_name <- as.character( deTab$chromosome_name );
  deTab$chromosome_name[ is.na( deTab$chromosome_name ) ] <- "NA";
  dDE <- deTab %>%
    filter( DE %in% c( -1, 1 ) ) %>% group_by( Contrast, chromosome_name ) %>% summarize( deNum = n() );
  dAll <- deTab %>% group_by( Contrast, chromosome_name ) %>% summarize( allNum = n() );
  d <- merge( dDE, dAll );
  ggplot( d, aes( x = allNum, y = deNum, label = chromosome_name, color = chromosome_name ) ) +
    geom_text( angle = -45 ) +
    guides( colour = FALSE ) +
    xlab( "Number of all genes" ) +
    ylab( "Number of DE genes" ) +
    facet_wrap( ~ Contrast, scales = "free_y" ) +
    theme_bw();
}

# ------------------------------------------------------------------------

volcanoPlot <- function( x,
  maxLabels = 50, labelVar = getTagVar( x ), colorVar = NULL,
  xVar = getDiffCountVar( x ),
  deAvgLine = NULL, allAvgLine = "smooth",
  linesAtFCs = c( 2 ), ...
) {
  deTab <- deData( x, ... );
  if( nrow( deTab ) == 0 ) { return( NULL ); }

  addManualColorScale <- FALSE;
  if( is.null( colorVar ) ) {
    deTab$color <- "gray";
    deTab$color[ deTab$DE != 0 & deTab$Rank > maxLabels ] <- "red";
    deTab$color[ deTab$DE != 0 & deTab$Rank <= maxLabels ] <- "black";
    colorVar <- "color";
    addManualColorScale <- TRUE;
  }

  deTab$minusLog10AdjPVal <- -log10( deTab$adj.P.Val );
  p <- ggplot( deTab, aes_string( x = xVar, y = "minusLog10AdjPVal", label = labelVar, color = colorVar ) ) +
    geom_hline( yintercept = -log10( 0.05 ), color = "blue", linetype = 2, alpha = 0.5 );
  for( fc in linesAtFCs ) {
    p <- p +
      geom_vline( xintercept = log2( fc ), color = "blue", linetype = 2, alpha = 0.5 ) +
      geom_vline( xintercept = -log2( fc ), color = "blue", linetype = 2, alpha = 0.5 );
  }

  ds <- subset( deTab, DE == 0 );
  if( nrow( ds ) > 0 ) p <- p +
    geom_point( data = ds, size = 1, alpha = 0.4 );

  p <- p +
    geom_vline( xintercept = 0, color = "darkgray", linetype = 2, alpha = 0.5 );

  ds <- subset( deTab, DE != 0 & ( Rank > maxLabels ) );
  if( nrow( ds ) > 0 ) p <- p +
    geom_point( data = ds, size = 1, alpha = 0.4 );

  ds <- subset( deTab, DE != 0 & ( Rank <= maxLabels ) );
  if( nrow( ds ) > 0 ) p <- p +
    geom_text( data = ds, size = 3, alpha = 1 );

  p <- p +
    theme_bw() +
    facet_wrap( ~ Contrast );

  if( addManualColorScale ) p <- p +
    scale_color_identity();

  p
}

# ------------------------------------------------------------------------

diffExpPlot <- function( x,
  maxLabels = 50, labelVar = getTagVar( x ),
  colorVar = NULL, selectVar = NULL,
  avgCountVar = getAvgCountVar( x ), diffCountVar = getDiffCountVar( x ),
  deAvgLine = NULL, allAvgLine = "smooth",
  linesAtFCs = c( 2 ),
  ...
) {
  deTab <- deData( x, ... );
  if( !is.null( selectVar ) ) {
    deTab <- deTab[ deTab[[ selectVar ]], ];
  }
  if( nrow( deTab ) == 0 ) { return( NULL ); }

  addManualColorScale <- FALSE;
  if( is.null( colorVar ) ) {
    deTab$color <- "gray";
    deTab$color[ deTab$DE != 0 & deTab$Rank > maxLabels ] <- "red";
    deTab$color[ deTab$DE != 0 & deTab$Rank <= maxLabels ] <- "black";
    colorVar <- "color";
    addManualColorScale <- TRUE;
  }

  p <- ggplot( deTab, aes_string( x = avgCountVar, y = diffCountVar, label = labelVar, color = colorVar ) );
  for( fc in linesAtFCs ) {
    p <- p +
      geom_hline( yintercept = log2( fc ), color = "blue", linetype = 2, alpha = 0.5 ) +
      geom_hline( yintercept = -log2( fc ), color = "blue", linetype = 2, alpha = 0.5 );
  }

  ds <- subset( deTab, DE == 0 );
  if( nrow( ds ) > 0 ) p <- p +
    geom_point( data = ds, size = 1, alpha = 0.4 );

  p <- p +
    geom_hline( yintercept = 0, color = "darkgray", linetype = 2, alpha = 0.5 );

  if( !is.null( allAvgLine ) ) {
    p <- p +
      stat_smooth( aes( group = 1 ), alpha = 0.3, fill = "lightblue", color = "lightblue" );
  }

  if( !is.null( deAvgLine ) ) {
    ds <- subset( deTab, DE != 0 );
    if( nrow( ds ) > 0 ) p <- p +
      stat_smooth( data = ds, aes( group = 1 ), alpha = 0.3, fill = "lightpink", color = "lightpink" );
  }

  ds <- subset( deTab, DE != 0 & ( Rank > maxLabels ) );
  if( nrow( ds ) > 0 ) p <- p +
    geom_point( data = ds, size = 1, alpha = 0.4 );

  ds <- subset( deTab, DE != 0 & ( Rank <= maxLabels ) );
  if( nrow( ds ) > 0 ) {
    if( "CI.L" %in% colnames( deTab ) ) p <- p +
        geom_errorbar( data = ds, aes( ymin = CI.L, ymax = CI.R ), alpha = 0.3 );
    p <- p +
      geom_text( data = ds, size = 3, alpha = 1 );
  }

  p <- p +
    theme_bw() +
    facet_wrap( ~ Contrast );

  if( addManualColorScale ) p <- p +
    scale_color_identity();

  p
}

# ------------------------------------------------------------------------

heatmapData <- function( x,
  maxAdjPVal = 0.05,
  countVar = getCountVar( x ),
  avgCountVar = getAvgCountVar( x ),
  diffVar = "diff",
  sampleVar = getSampleVar( x ),
  tagVar = getTagVar( x ),
  ...
) {
  de <- deData( x, ... );
  de <- de %>%
    group_by_( .dots = tagVar ) %>%
    filter( min( adj.P.Val ) <= maxAdjPVal ) %>%
    ungroup();

  dc <- countsData( x, addSampleSheet = TRUE ) %>%
    semi_join( de, by = tagVar );
  dc[[ diffVar ]] <- dc[[ countVar ]] - dc[[ avgCountVar ]];
  dc <- dc %>%
    select_( .dots = c( tagVar, sampleVar, diffVar ) ) %>%
    spread_( sampleVar, diffVar );
  #samples <- as.character( getSampleSheet( x )[[ sampleVar ]] );
  #m <- m[ , samples ];

  dc
}

heatmapPlot <- function( x, ..., diffVar = "diff", labRow = NA, Colv = NA ) {
  d <- heatmapData( x, diffVar = diffVar, ... );

  tagVar <- getTagVar( x );
  m <- d %>%
    select_( .dots = setdiff( colnames( d ), tagVar ) ) %>%
    as.matrix();
  rownames( m ) <- d[[ tagVar ]];

  heatmap3::heatmap3( m, labRow = labRow, balanceColor = TRUE, Colv = Colv );
}

# ------------------------------------------------------------------------

topTagsData <- function( x,
  topTagsNum = 20,  # obsoleting, use maxRank instead; will be removed after 20161001,
  contrastNames = colnames( getContrasts( x ) ),
  DEs = c( 1, -1 ),
  maxRank = topTagsNum,
  onlyContrastMatch = TRUE,
  ...
) {
  sampleVar <- getSampleVar( x );
  tagVar <- getTagVar( x );

  de <- deData( x, ... ) %>%
    filter( Rank <= maxRank ) %>%
    filter( Contrast %in% contrastNames ) %>%
    filter( DE %in% DEs );
  dc <- countsDensityData( x );
  d <- de %>% left_join( dc, by = tagVar )  %>% arrange( Contrast, Rank )

  if( onlyContrastMatch ) {
    m <- as.data.frame( getDesign( x ) %*% getContrasts( x ) ) %>% add_rownames( var = sampleVar );

    needRefContrast <- names( which( colSums( getContrasts( x ) ) != 0 ) );
    m <- do.call( rbind, lapply( needRefContrast, function( cn ) {
      m$idx <- m %>% group_indices_( .dots = setdiff( needRefContrast, cn ) );
      selIdxs <- ( m %>% select_( .dots = c( cn, "idx" ) ) %>% distinct() %>% group_by( idx ) %>%
        summarize( count = n() ) %>% filter( count == max( count ) ) )$idx;
      m <- m %>% filter( idx %in% selIdxs ) %>% select_( .dots = sampleVar );
      m$Contrast <- cn;
      m
    } ) )

    d %>% semi_join( m )
  } else {
    d
  }
}

topTagsPlot <- function( x,
  ...,
  coordVar = getTagVar( x ),
  countVar = getCountVar( x ),
  colorGroup = FALSE, # obsoleting, use colorVar instead; will be removed after 20161001
  colorVar = colorGroup,
  shapeVar = NULL,
  alpha = 1,
  dataFun = function( x, ... ){ topTagsData( x, ... ) }
) {
  d <- dataFun( x, ... );
  groups <- getGroups( x );
  if( is.logical( colorVar ) ) {
    colorVar <- getGroupVar( x );
    if( is.null( colorVar ) ) colorVar <- .getGroupName( x, 1 );
  }

  d$ContrastRank <- paste0( d$Contrast, ":", d$Rank );
  dd <- d %>% dplyr::select( Contrast, Rank, ContrastRank ) %>% distinct() %>% arrange( Contrast, desc( Rank ) );
  d$ContrastRank <- factor( d$ContrastRank, levels = as.character( dd$ContrastRank ) );

  if( nrow( d ) > 0 ) {
    p <- ggplot( d, aes_string(
      x = countVar,
      y = "ContrastRank",
      color = if( !is.null( colorVar ) ) colorVar else NULL,
      shape = if( !is.null( shapeVar ) ) shapeVar else NULL
    ) );
    if( !is.null( shapeVar ) ) {
      p <- p + geom_point( size = 2, alpha = alpha );
    } else {
      p <- p + geom_point( size = 5, alpha = alpha, shape = "|" );
    }
    p <- p +
      facet_grid( Contrast ~ ., scales = "free_y" ) +
      theme_bw();
    if( !is.null( colorVar ) ) {
      colorValues <- groups[[ colorVar ]]$df$color;
      names( colorValues ) <- groups[[ colorVar ]]$df$level;
      p <- p + scale_color_manual( name = groups[[ colorVar ]]$name, values = colorValues );
    }
    if( coordVar != "ContrastRank" ) {
      p <- p + scale_y_discrete( breaks = d$ContrastRank, labels = d[[ coordVar ]], name = coordVar );
    }
    p
  } else {
    NULL
  }
}

# ------------------------------------------------------------------------

# topTags.DETable->deData: consider which arguments are needed maxRank = 50, topTagsNum = maxRank, contrastNames = colnames( getContrasts( x ) ), DEs = c( 1, -1 )
# will be replaced by pairedTopTagsPlot
plotPairedTopTags <- function( deTable, ...,
                         subjectVar = getSubjectVar( deTable ),
                         groupVar = getGroupVar( deTable ),
                         groupLevels = as.character( head( ( getGroups( deTable ) )[[ groupVar ]]$df$level, 2 ) ),
                         countName = "log2CPM"
) {
  d <- topTags.DETable( deTable, ... );
  counts <- getCounts( deTable );
  groups <- getGroups( deTable );

  df <- melt( counts, value.name = countName );
  if( is.numeric( df$sample ) ) df$sample <- factor( df$sample );
  tagName <- colnames( df )[[1]];
  tagVar <- getTagVar( deTable );
#  df <- merge( df, d[ c( "tag", "Rank", "DE", "Contrast" ) ], by.x = tagName, by.y = "tag", all.x = FALSE, all.y = TRUE );
  df <- merge( d[ c( tagVar, "Rank", "DE", "Contrast" ) ], df, by.y = tagName, by.x = tagVar, all.y = FALSE, all.x = TRUE );
  df$Rank <- factor( df$Rank, levels = sort( unique( df$Rank ), decreasing = FALSE ) );
  df[[ groupVar ]] <- groups[[ groupVar ]]$sf[ df$sample, ]$level;
  df[[ subjectVar ]] <- groups[[ subjectVar ]]$sf[ df$sample, ]$level;
  dx <- df[ df[[ groupVar ]] %in% groupLevels[[1]], ];
  dy <- df[ df[[ groupVar ]] %in% groupLevels[[2]], ];
  rn <- range( dx[[ countName ]], dy[[ countName ]], na.rm = TRUE );
  df <- merge( dx, dy, by = setdiff( colnames( dx ), c( "sample", countName, groupVar ) ) );

  if( !is.null( subjectVar ) ) {
    df$color <- groups[[ subjectVar ]]$df[ df[[ subjectVar ]], ]$level;
  }
  if( nrow( df ) > 0 ) {
    p <- ggplot( df, aes_string(
      x = paste0( countName, ".x" ),
      y = paste0( countName, ".y" ),
      color = if( !is.null( subjectVar ) ) "color" else NULL
    ) ) +
      geom_abline( slope = 1, color = "gray", linetype = 3 ) +
      geom_abline( slope = 1, intercept = +1, color = "blue", linetype = 2, alpha = 0.5 ) +
      geom_abline( slope = 1, intercept = -1, color = "blue", linetype = 2, alpha = 0.5 ) +
      geom_point( alpha = 0.7 ) +
      theme_bw() + xlim( rn ) + ylim( rn ) +
      xlab( paste0( countName, " (", groupVar, " == ", groupLevels[[1]], ")" )) +
      ylab( paste0( countName, " (", groupVar, " == ", groupLevels[[2]], ")" ));
    if( nlevels( df$Contrast ) > 1 ) {
      p <- p + facet_wrap( ~ Contrast + Rank + tag );
    } else {
      p <- p + facet_wrap( ~ Rank + tag );
    }
    if( !is.null( subjectVar ) ) {
      colorValues <- groups[[ subjectVar ]]$df$color;
      names( colorValues ) <- groups[[ subjectVar ]]$df$level;
      p <- p + scale_color_manual( name = groups[[ subjectVar ]]$name, values = colorValues );
    }
    p
  } else {
    NULL
  }
}

# ------------------------------------------------------------------------

if( 0 ) { # obsoleted 20160613
  plotUpDownBalance <- function( deTable, ...,
                                 topTagsNum = +Inf,
                                 coordVar = "Rank",
                                 decreasing = FALSE
  ) {
    if( coordVar == "logFC" ) coordVar <- "log2FC";
    if( coordVar == "AveExpr" ) coordVar <- "log2AvgExpr";

    d <- topTags.DETable( deTable, topTagsNum = topTagsNum, ... );
    if( nrow( d ) > 0 ) {
      d <- d[ order( d[[ coordVar ]], decreasing = decreasing ), ];
      d <- do.call( rbind, lapply( split( d, d$Contrast ), function( d ) {
        d$Upregulated <- cumsum( d$DE == 1 );
        d$Downregulated <- cumsum( d$DE == -1 );
        d
      } ) );
      div <- round( nrow( d ) / 10 ) + 1;
      dd <- melt(
        d[ c( "log2AvgExpr", "Contrast", "Rank", "Upregulated", "Downregulated" ) ],
        id.vars = c( "log2AvgExpr", "Contrast", "Rank" ),
        variable.name = "Direction",
        value.name = "Count"
      );
      p <- ggplot( dd, aes_string(
        x = coordVar, y = "Count", group = "Direction", color = "Direction", shape = "Direction"
      )) +
        geom_line( alpha = 0.5, linetype = 1 ) +
        geom_point( data = subset( dd, ( as.integer( Rank ) %% div ) == 0 ), size = 4 ) +
        theme_bw() +
        facet_grid( ~ Contrast ) +
        scale_shape_manual( values = c( Upregulated = 2, Downregulated = 6 ) )
      p
    } else {
      NULL
    }
  }
}

## -----------------------------------------------------------------------
## deComparisonPlot
## -----------------------------------------------------------------------

deComparisonTable <- function( deTable.x, deTable.y, ... ) {
  ret <- merge( deData( deTable.x ), deData( deTable.y ), by = c( "tag", "Contrast" ), ... );
  class( ret ) <- c( "DEComparisonTable", class( ret ) );
  ret
}

plot.DEComparisonTable <- function( x,
                                    coordVar = "log2FC", maxLabels = 30, labelVar = "tag",
                                    maxColorPoints = maxLabels, deFacetGrid = TRUE,
                                    contrastNames = NULL, xLab = "x", yLab = "y"
) {
  if( !is.null( contrastNames ) ) x <- subset( x, Contrast %in% contrastNames );

  p <- ggplot( x, aes_string(
    x = paste0( coordVar, ".x" ),
    y = paste0( coordVar, ".y" ),
    label = labelVar
  ) ) +
    #    geom_bin2d( alpha = 0.5 ) +
    geom_abline( slope = 1, color = "gray", linetype = 2 ) +
    geom_point( data = subset( x, !( Rank.x <= maxColorPoints & DE.x != 0 ) & !( Rank.y <= maxColorPoints & DE.y != 0 ) ), alpha = 0.3, color = "gray", shape = 20, size = 1 );
  p <- p +
    geom_point( data = subset( x, ( Rank.x <= maxColorPoints & DE.x != 0 ) & !( Rank.y <= maxColorPoints & DE.y != 0 ) ), alpha = 0.6, color = "blue", size = 2 ) +
    geom_point( data = subset( x, !( Rank.x <= maxColorPoints & DE.x != 0 ) & ( Rank.y <= maxColorPoints & DE.y != 0 ) ), alpha = 0.8, color = "violet", size = 2 ) +
    geom_point( data = subset( x, ( Rank.x <= maxColorPoints & DE.x != 0 ) & !( Rank.y <= maxColorPoints & DE.y != 0 ) ), alpha = 1, color = "blue", size = 4, shape = "|" ) +
    geom_point( data = subset( x, !( Rank.x <= maxColorPoints & DE.x != 0 ) & ( Rank.y <= maxColorPoints & DE.y != 0 ) ), alpha = 1, color = "violet", size = 13, shape = "-" ) +
    geom_point( data = subset( x, ( Rank.x <= maxColorPoints & DE.x != 0 ) & ( Rank.y <= maxColorPoints & DE.y != 0 ) ), alpha = 0.8, color = "red", size = 2 ) +
    geom_text( data = subset( x, ( Rank.x <= maxLabels & DE.x != 0 ) & !( Rank.y <= maxLabels & DE.y != 0 ) ), alpha = 1, color = "blue", size = 3 ) +
    geom_text( data = subset( x, !( Rank.x <= maxLabels & DE.x != 0 ) & ( Rank.y <= maxLabels & DE.y != 0 ) ), alpha = 1, color = "violet", size = 3 ) +
  geom_text( data = subset( x, ( Rank.x <= maxLabels & DE.x != 0 ) & ( Rank.y <= maxLabels & DE.y != 0 ) ), alpha = 1, color = "red", size = 3 );
    p <- p +
    xlab( paste0( coordVar, " [", xLab, "]" ) ) +
    ylab( paste0( coordVar, " [", yLab, "]" ) ) +
    theme_bw();
  if( deFacetGrid ) p <- p + facet_grid( DE.y ~ DE.x, labeller = label_both );
  p
}

deComparisonPlot <- function( deTable.x, deTable.y, ... ) {
  m <- deComparisonTable( deTable.x, deTable.y );
  plot( m, ... );
}

## -----------------------------------------------------------------------
## DEComparison
## -----------------------------------------------------------------------

# byContrasts = length( l ) == 1, byDEs = length( l ) > 1
compareDEs <- function( ... ) {
  l <- list( ... );
  if( is.null( names( l ) ) ) names( l ) <- paste0( "DE", 1:length( l ) );
  cols <- c( "tag", "DE", "Contrast" );
  d <- do.call( rbind, lapply( names( l ), function( n ) {
    data.frame( deData( l[[n]] )[ cols ], method = n )
  } ) );

  dd <- dcast(
    d[ c( "Contrast", "tag", "method", "DE" ) ],
    Contrast + tag ~ method,
    value.var = "DE"
  );

  deByMethod <- d %>% group_by( Contrast, method, DE ) %>%
    summarize( count = n() ) %>% ungroup() %>%
    dcast( Contrast + method ~ DE, value.var = "count", fill = 0 );

  dd <- d %>% dcast( tag ~ method + Contrast, value.var = "DE" );
  nm <- colnames( dd %>% dplyr::select( -tag ) );
  deByConsistency <- dd %>% group_by_( .dots = nm ) %>%
    summarise( count = n() ) %>% ungroup() %>% arrange( desc( count ) );

  if( 0 ) {
    combn( nm, 2, function( nm ) {
      dd <- dd %>% select_( .dots = nm );
      dd %>% group_by_( .dots = colnames( dd ) ) %>% summarise( count = n() ) %>%
        ungroup() %>% arrange( desc( count ) )
    }, simplify = FALSE )
  }

  ret <- list(
    deList = l,
    deDEs = dd,
    deByMethod = deByMethod,
    deByConsistency = deByConsistency
  )
  class( ret ) <- c( "DEComparison", class( ret ) );
  ret
}

summary.DEComparison <- function( x ) {
  x[ c( "deByMethod", "deByConsistency" ) ];
}

writeReport.DEComparison <- function( x, names = c(), ... ) {
  for( n in names( x$deList ) ) {
    writeReport( x$deList[[ n ]], names = c( names, n ), ... );
  }

  writeReport( x$deByMethod, names = c( names, report = "deByMethod" ), ... );

  writeReport( x$deDEs, names = c( names, report = "deDEs" ), ... );

  if( 0 ) {
    deDEs <- x$deDEs;
    for( d in split( deDEs, factor( deDEs$Contrast ) ) ) {
      if( nrow( d ) > 0 ) {
        contrastName <- as.character( d$Contrast[[1]] );
        writeReport( d, names = c( names, report = "deByGenes", contrast = contrastName ), ... );
      }
    }
  }

  writeReport( x$deByConsistency, names = c( names, report = "deByConsistency" ), ... );
}

plot.DEComparison <- function( x, methods = c( 1, 2 ), ... ) {
  ll <- x$deList[ methods ];
  deComparisonPlot( ll[[1]], ll[[2]], xLab = names( ll )[[1]], yLab = names( ll )[[2]], ... );
}

vennPlot.DEComparison <- function( x, DE, contrastNameSuffix = NULL ) {
  d <- x$deDEs %>% dplyr::select_( .dots = colnames( x$deDEs )[-1] ) %>% as.matrix;
  mode( d ) <- "numeric";
  d[ is.na( d ) ] <- 0;
  if( !is.null( contrastNameSuffix ) ) {
    d <- d[ , grepl( contrastNameSuffix, colnames( d ) ) ];
  }

  stopifnot( !any( apply( d, 1, function( v ) any( v > 0 ) && any( v < 0 ) ) ) );

#  library( venneuler );
#  plot( venneuler( d[ noMinusRows, ] != 0 ) )

#  vennDiagram( vennCounts(d[ noMinusRows, ] != 0 ) )

  library( gplots );

  if( DE == +1 ) {
    noMinusRows <- apply( d, 1, function( v ) all( v >= 0 ) && any( v > 0 ) );
    venn( as.data.frame( d[ noMinusRows, ] != 0 ) );
  } else if( DE == -1 ) {
    noPlusRows <- apply( d, 1, function( v ) all( v <= 0 ) && any( v < 0 ) );
    venn( as.data.frame( d[ noPlusRows, ] != 0 ) );
  } else {
    stop( "DE must be either +1 or -1." )
  }
}

## -----------------------------------------------------------------------

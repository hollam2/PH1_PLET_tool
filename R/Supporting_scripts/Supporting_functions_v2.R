#functions for PH1/FW5 code

#function for remove years from time series with less than n months of interpolated data and determine proportion of years removed
clean_years <- function(x, thr){
  temp <- x %>%
    group_by(lifeform) %>%
    dplyr::mutate(orig_year = length(unique(year))) %>% #count number of years per cell in the original data
    ungroup() %>%
    filter(!is.na(abundance)) %>%
    group_by(lifeform, year) %>% 
    dplyr::mutate(n=n()) %>% #count number of months of data within each year per cell
    ungroup() %>%
    filter(n < thr) %>% #extract years with less than n months of real data
    group_by(lifeform) %>%
    dplyr::mutate(prop_years_removed = length(unique(year)) / orig_year) %>% #determine proportion of years of data removed per grid cell
    ungroup()
  
  if(nrow(temp > 0)){
    print(paste("Proportion of years removed:", unique(round(temp$prop_years_removed, 2))))
  } else {
    print(paste("Proportion of years removed:", 0))
  }
  
  years_exclude <- sort(unique(temp$year))
  
  output <- x
  output[output$year %in% years_exclude,"abundance"] <- NA

  return(output)
}

#function for filling month gaps in the time series
fill_gaps <- function(x, max_gap = 3){
  
  #arrange the data to allow for NA interpolation across time
  temp <- x %>% 
    arrange(lifeform, year, month) %>%
    group_by(lifeform) %>%
    dplyr::mutate(abundance = zoo::na.approx(abundance, maxgap = max_gap, rule = 2)) %>%
    ungroup()
  
  print(paste("Proportion of missing months filled:", round(1-(sum(is.na(temp$abundance))/sum(is.na(x$abundance))),2)))
  
  output <- temp %>%
    filter(!is.na(abundance))
  
  return(output)
}

#construct a dataframe of relevant lifeform pair comparisons
df_lf <- rbind(data.frame(V1 = "diatom", V2 = "dinoflagellate"),
               data.frame(V1 = "tycho_diatoms", V2 = "pelagic_diatoms"),
               data.frame(V1 = "lg_copepods", V2 = "sm_copepods"),
               data.frame(V1 = "holoplankton", V2 = "meroplankton"),
               data.frame(V1 = "lg_phyto", V2 = "sm_phyto"),
               data.frame(V1 = "phytoplankton", V2 = "noncarniv"),
               data.frame(V1 = "crustacean", V2 = "gelatinous"),
               data.frame(V1 = "gelatinous", V2 = "fishlarvae")
)

#function for extracting a dataframe for a particular time period
dataSelect <- function(x, lf, lims){
  
  output <- x %>%
    filter(lifeform %in% as.vector(unlist(lf)),
           year>=lims[1] & year<=lims[2]) %>%
    ungroup() %>%
    pivot_wider(names_from = lifeform, values_from = abundance)
  
  return(output)
}


#function to prepare the reference envelopes for the multiple lifeform pairs comparisons
find_envAll <- function(x, lf){
  
  #determine relevant lifeform pairs for the dataset
  x_temp <- x
  
  #find the relevant lifeform pairs
  x_temp <- x_temp[,colSums(is.na(x_temp))<nrow(x_temp)]
  
  lf_temp <- lf %>%
    filter(V1 %in% all_of(colnames(x_temp)),
           V2 %in% all_of(colnames(x_temp)))
  
  main_outer <- data.frame()
  main_inner <- data.frame()
  for(i in 1:nrow(lf_temp)){
    temp_lf <- as.vector(unlist(lf_temp[i,]))
    
    temp_x <- x_temp %>%
      dplyr::select(all_of(temp_lf))

    df_outer <- data.frame()
    df_inner <- data.frame()
    
      #command to skip envelope fitting for data with no variance
      abort <- ifelse(sd(as.vector(unlist(temp_x[,1]))) == 0 | sd(as.vector(unlist(temp_x[,2])))==0, TRUE, FALSE)
      
      if(abort==FALSE){
        
        envPts <- findEvn(as.vector(unlist(temp_x[,1])),
                          as.vector(unlist(temp_x[,2])),
                          p=0.9,
                          sc=TRUE)
        envPts_unlist <- rbindlist(envPts, fill=TRUE)
        temp_outer <- data.frame(outX=envPts_unlist$outX[complete.cases(envPts_unlist$outX)],
                                 outY=envPts_unlist$outY[complete.cases(envPts_unlist$outY)])
        temp_inner <- data.frame(inX=envPts_unlist$inX[complete.cases(envPts_unlist$inX)],
                                 inY=envPts_unlist$inY[complete.cases(envPts_unlist$inY)])
        
        df_outer <- rbind(df_outer, temp_outer)
        df_inner <- rbind(df_inner, temp_inner)
      }
    
    
    if(nrow(df_outer) > 0 & nrow(df_inner) > 0){
      
      df_outer$lf_pair <- paste(temp_lf, collapse="-")
      df_inner$lf_pair <- paste(temp_lf, collapse="-")
      
      main_outer <- rbind(main_outer, df_outer)
      main_inner <- rbind(main_inner, df_inner)
      
    
    }
}
  
  main_list <- list(main_outer, main_inner)
  return(main_list)
}


#function to calculate the lifeform pairs indicator from the reference envelopes and comparison data
PIcalcAll <- function(x, y, z, lf){
  
  main_outer <- x[[1]]
  main_inner <- x[[2]]
  
  main_output <- data.frame()
  for(i in 1:length(unique(main_outer$lf_pair))){
    temp_lf <- unlist(strsplit(sort(unique(main_outer$lf_pair))[i], "-"))
    
    temp_outer <- subset(main_outer, lf_pair == paste(temp_lf, collapse="-"))
    temp_inner <- subset(main_inner, lf_pair == paste(temp_lf, collapse="-"))
    
    temp_y <- y %>%
      dplyr::select(all_of(temp_lf))
    
    df_z <- z %>%
      dplyr::select(all_of(temp_lf))

    #arrange the envelope data back into a list
    envelopePts <- list("EnvOuter"=data.frame("outX" = temp_outer$outX,"outY"=temp_outer$outY),
                          "EnvInner"=data.frame("inX" = temp_inner$inX,"inY" = temp_inner$inY))
      
    compDat <- data.frame(y1 = as.vector(unlist(temp_y[,1])),
                            y2 = as.vector(unlist(temp_y[,2])))
      
    #add labelling variables to PI results dataframe
    df_refPoints <- df_z %>%
        dplyr::summarise(refPoints = n())
      
    #command to skip envelope fitting for data with no reference envelope
    abort <- ifelse(nrow(envelopePts[[1]]) == 0 & nrow(envelopePts[[2]]) == 0, TRUE, FALSE)
      
      if(abort==FALSE){
        piPts <- PIcalc(compDat, envelopePts, 0.9)
        piPts <- do.call(cbind.data.frame, piPts)
        piPts$refPoints <- df_refPoints$refPoints[1]
        piPts$lf_pair <- paste(temp_lf, collapse="-")
      }
    
    main_output <- rbind(main_output, piPts)
  }
  
  main_output <- dplyr::rename(main_output, binomialProbability = 'binomial probability')
  
  return(main_output)
}


#function for plotting the PI envelope
plot_env <- function(x, y, z, lf, pi){
  
  #rounding function for labelling
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
  
  refs <- y %>%
    dplyr::summarise(refStart=min(year),
                     refStop=max(year))
  
  comps <- z %>%
    dplyr::summarise(compStart=min(year),
                     compStop=max(year))

  
  df_lookup_main <- data.frame(lf_pair=pi$lf_pair) %>%
   mutate(refStart=refs$refStart[1],
          refStop=refs$refStop[1],
          compStart=comps$compStart[1],
          compStop=comps$compStop[1]) %>%
    dplyr::mutate(years_label = paste0("Ref: " , refStart, "-", refStop, ",", " Comp: ", compStart, "-", compStop))
  
  #create labeller lookup table
  df_lookup_main$string <- paste0('PI: ', specify_decimal(pi$PI, 2), ' ',
                                  df_lookup_main$years_label,'\n',
                                  'ref points: ', pi$refPoints, ',', ' ',
                                  'comp points: ', pi$newPoints,',', ' ',
                                  'binom-p: ', specify_decimal(pi$binomialProbability, 4), ',', ' ',
                                  'chi-sq: ', specify_decimal(pi$chi.sq, 1))
  
  main_outer <- x[[1]]
  main_inner <- x[[2]]
  
  plot_list <- list()
  for(i in 1:length(unique(main_outer$lf_pair))){
    
    temp_lf <- unlist(strsplit(sort(unique(main_outer$lf_pair))[i], "-"))
    
    df_outer <- subset(main_outer, lf_pair == paste(temp_lf, collapse="-"))
    df_inner <- subset(main_inner, lf_pair == paste(temp_lf, collapse="-"))
    
    names(df_outer)[1:2] <- c("x", "y")
    names(df_inner)[1:2] <- c("x", "y")
    
    df_outer$subid <- 1L
    df_inner$subid <- 2L
    df_polys <- rbind(df_outer, df_inner)
    
    temp_ref <- y %>% dplyr::select(1:4, temp_lf) %>%
      arrange(year, month)
    names(temp_ref)[c(ncol(temp_ref)-1,ncol(temp_ref))] <- c("vx", "vy")
    
    temp_comp <- z %>% dplyr::select(1:4, temp_lf) %>%
      arrange(year, month)
    names(temp_comp)[c(ncol(temp_comp)-1,ncol(temp_comp))] <- c("vx", "vy")
    
    #check the proportion of years represented in the data
    year <- rbind(temp_comp, temp_ref) %>%
      dplyr::select(year) %>%
      distinct() %>%
      dplyr::mutate(prop_years = length(unique(year) %in% seq(min(year), max(year), 1)) / length(seq(min(year), max(year), 1)))
    
    #grouping factor for colouring months
    temp_comp$month <- as.numeric(temp_comp$month)
    temp_comp$season <- ifelse(temp_comp$month %in% c(1, 2, 12), "months: 1 2 12",
                               ifelse(temp_comp$month %in% c(3, 4, 5), "months: 3 4 5",
                                      ifelse(temp_comp$month %in% c(6, 7, 8), "months: 6 7 8",
                                             ifelse(temp_comp$month %in% c(9, 10, 11), "months: 9 10 11", "ERROR"))))
    
    #Create a custom color scale
    factor_levels <- unique(temp_comp$season)
    myColors <- c("blue", "green", "yellow", "red")
    names(myColors) <- levels(factor_levels)
    
    #subset lookup table to panel of relevance
    df_lookup_temp <- unique(subset(df_lookup_main, lf_pair == paste(temp_lf, collapse="-")))
    
    #create title reference string
    years_label <- paste0("Ref: " , df_lookup_temp$refStart[1], "-", df_lookup_temp$refStop[1], ",", 
                          " Comp: ", df_lookup_temp$compStart[1], "-", df_lookup_temp$compStop[1])
    df_lookup_temp <- setNames(df_lookup_temp$string, df_lookup_temp$lf_pair)
    
      gg_panel <- ggplot() +
        geom_polygon(data=df_polys, aes(x, y, subgroup = subid), fill="grey60", colour="black", size=0.25, alpha=0.5) +
        geom_path(data=temp_comp, aes(x=vx, y=vy), colour="grey", linetype = 2, size=0.25) +
        geom_point(data=temp_comp, aes(x=vx, y=vy, fill=season), shape=21) +
        geom_point(data=temp_ref, aes(x=vx, y=vy), shape=21, fill=NA, colour="grey60", alpha=0.5) +
        scale_x_continuous(name=bquote(log[10]* "(" * .(temp_lf[1]) * ")")) +
        scale_y_continuous(name=bquote(log[10]* "(" * .(temp_lf[2]) * ")")) +
        scale_fill_manual(values=myColors)+
        facet_wrap(~ lf_pair, scales="free", labeller = labeller(lf_pair=df_lookup_temp)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              legend.title = element_blank(),
              legend.position = "none")
      
    plot_list[[paste(temp_lf, collapse="-")]] <- gg_panel
    
  }
  return(plot_list)
}


#function for modelling change in lifeforms over time with Kendall test
kendallAll <- function(x){
  
  #generate relevant dataset within year limits
  trajData <- x
  
  #model annual change in abundance of each group in each spatial unit
  df_fits_tot <- trajData %>%
    dplyr::select(year, month, num_samples, lifeform, abundance) %>%
    group_by(year, lifeform) %>%
    dplyr::summarise(abundance_mean = mean(abundance, na.rm=T),
                     .groups="drop") %>%
    filter(!is.nan(abundance_mean)) %>%
    group_by(lifeform) %>%
    dplyr::mutate(count = n()) %>%
    filter(count >= 3) %>%
    dplyr::mutate(year=as.integer(year)) %>%
    nest() %>%
    mutate(fits = map(data, ~EnvStats::kendallTrendTest(abundance_mean ~ year, ci.slope=FALSE, data=.x)),
           fits2 = map(fits, ~structure(.x, class="htest")),
           fits3 = map(fits2, ~tidy(.x))) %>%
    dplyr::select(-c(data, fits, fits2)) %>%
    unnest_wider(fits3)
  
  #column to code whether p-value is significant
  df_fits_tot$sig <- ifelse(df_fits_tot$p.value <= 0.05, TRUE, FALSE)
  
  #simplify the dataframe output
  df_fits_tot <- data.frame(lifeform = df_fits_tot$lifeform,
                            statistic = df_fits_tot$statistic,
                            p = df_fits_tot$p.value,
                            sig = df_fits_tot$sig) %>%
    arrange(lifeform)
  
  
  
  return(df_fits_tot)
}


#function to prepare the data to be plotted as time-series
create_ts <- function(x, y, limits, ids){
  
  #generate relevant dataset within year limits
  trajData <- x
  
  #match model results to new plotting dataframe
  df_plot <- merge(trajData, y, by=c("lifeform"))
  
  #convert year and month variables to date format
  df_plot$date_mon <- as.Date(paste(df_plot$year, df_plot$month, 16), "%Y %m %d")
  df_plot$date_year <- as.Date(paste(df_plot$year, 07, 02), "%Y %m %d")
  
  #tally number of samples for each facet
  df_plot <- df_plot %>%
    dplyr::mutate(year=as.integer(year)) %>%
    filter(!is.na(abundance)) %>%
    group_by(lifeform) %>%
    dplyr::mutate(sumSamples = sum(num_samples, na.rm=T),
                  prop_years = length(unique(year) %in% seq(min(year), max(year), 1)) / length(seq(min(year), max(year), 1))) %>%
    ungroup()

  return(df_plot)
}


#function for plotting time-series
plot_ts <- function(x){
  
  #create labeller lookup table
  df_lookup_main <- data.frame(lifeform=x$lifeform) %>%
    dplyr::mutate(string = paste0('z: ', round(x$statistic, 2), '   ', 
                                  'p: ', ifelse(x$p <= 0.05, "<=0.05", round(x$p,3)), '  ',
                                  'n: ', x$sumSamples)) %>%
    distinct()
  
  
  #generate a temp dataframe
  x_temp <- x
  
  plot_list <- list()
  for (i in 1:length(unique(x_temp$lifeform))){
    
    lf_temp <- sort(unique(x_temp$lifeform))[i]
    
    #subset to polygon of interest
    df_lookup_temp <- unique(subset(df_lookup_main, lifeform == lf_temp))
    
    #filter to polygon of interest and plot results 
    temp <- x_temp %>%
      filter(lifeform == lf_temp)
    
    df_lookup_temp_id <- setNames(df_lookup_temp$string, df_lookup_temp$lifeform)
      
      temp_params <- temp %>%
        dplyr::select(year, abundance)
      
      years <- c(min(temp_params$year),max(temp_params$year))
      years_brk <- ifelse(plyr::round_any(years[2],5,f=ceiling)-plyr::round_any(years[1],5,f=floor) <= 20, "2 years", "5 years")
      
      years <- seq.Date(from = as.Date(paste(plyr::round_any(years[1],5,f=floor),"01","01",sep="-")), 
                        to = as.Date(paste(plyr::round_any(years[2],5,f=ceiling)+1,"01","01",sep="-")), 
                        by = years_brk)
      
      y_bks <- if(max(temp_params$abundance)>=3){
        seq(0,10,1)
      }else if(max(temp_params$abundance)<3 & max(temp_params$abundance)>1){
        seq(0,10,0.5)  
      }else{
        seq(0,10,0.2)
      }
      
      gg_panel <- temp %>%
        group_by(year) %>%
        dplyr::mutate(abundance_annual = mean(abundance, na.rm=T)) %>%
        ungroup() %>%
        dplyr::mutate(interp=ifelse(num_samples==0,TRUE,FALSE)) %>%
        arrange(year, month) %>%
        ggplot(.,aes(date_mon, abundance, colour=interp))+
        geom_path(aes(group=1),size=0.25)+
        geom_smooth(aes(date_year, abundance_annual), formula= y ~ x,
                    linetype="dashed", colour="black", 
                    method = 'lm', se = FALSE) +
        geom_line(aes(date_year, abundance_annual), colour="blue", size=1)+
        geom_point(aes(date_year, abundance_annual), shape=21, fill="blue")+
        facet_wrap(~lifeform, ncol=1, scales="free_y", labeller = labeller(lifeform=df_lookup_temp_id))+
        scale_x_date(breaks = years, date_labels = "%Y",
                     minor_breaks = NULL)+
        scale_y_continuous(name=bquote(log[10]* "(" * .(lf_temp) * ")"), minor_breaks = NULL, breaks=y_bks)+
        scale_colour_manual(values=c("TRUE"="grey","FALSE"="blue"))+
        guides(shape="none")+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
              axis.title.x = element_blank(),
              legend.position = "none")
      
    plot_list[[lf_temp]] <- gg_panel
  }
  return(plot_list)
}


#function to select, combine and save the combined plots
combine_pi_plots <- function(x, y, z, limits, path){
  
  #generate directory for the plots
  output_path_traj <- paste(path, "timeseries_", limits[1], "_to_", limits[2], "/", sep="")
  dir.create(file.path(output_path_traj), showWarnings = FALSE)
  do.call(file.remove, list(list.files(output_path_traj, full.names = TRUE)))
  
  for(i in 1:length(x)){
    
    lf1 <- unlist(strsplit(names(x)[i], "-"))[1]
    lf2 <- unlist(strsplit(names(x)[i], "-"))[2]
    
    lf_pair_temp <- paste0(lf1,"-", lf2) 
    lf_pair_plot <- x[[lf_pair_temp]]
    
    lf1_plot <- y[[lf1]]
    lf2_plot <- y[[lf2]]
    
    #print an output to provide the user with a sense of progress
    print(paste0("lifeform pair: ", paste0(lf1,"-", lf2)))
      
    temp_plot <- grid.arrange(lf_pair_plot, grid.arrange(lf1_plot, lf2_plot), 
                                nrow = 1, 
                                widths=c(1,2.5))
      
    #create the filename
    filename_temp <- paste0(output_path_traj, paste0(lf1,"-", lf2,".png"))
      
    ggsave(temp_plot, file=filename_temp,
            height=15, width=45, units="cm", bg="white")

  }
}

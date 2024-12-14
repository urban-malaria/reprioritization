source("functions.R", echo=F)




shapefile <- lapply(seq_along(reprioritization_states), 
                    function(rowz) sf::st_read(file.path(shapefile, "STATES", reprioritization_states[rowz],
                                                               paste0(reprioritization_states[rowz],"_State.shp"))))



ward_data <- lapply(seq_along(reprioritization_states), 
                    function(rowz) read.csv(file.path(datafile,"five_extractions", 
                                                      paste0(reprioritization_states[rowz],
                                                             "_wards_variables.csv"))))



variable_relationships <- data.frame(mean_EVI = "direct", 
                                     mean_NDVI = "direct",
                                     mean_rainfall = "direct", 
                                     distance_to_water = "inverse", 
                                     RH_mean = "direct", 
                                     temp_mean  = "direct",
                                     housing_quality = "direct", 
                                     pfpr = "direct", 
                                     tpr = "direct", 
                                     avgRAD = "direct", 
                                     settlement_type = "direct"
                                     ,
                                     flood = "direct",
                                     NDWI = "direct",
                                     NDMI  = "direct"
                                     )


possible_vars <- c("mean_EVI", "housing_quality", "distance_to_water",
                   "settlement_type","pfpr", "tpr"
                   ,
                   "flood",
                   "NDWI", "NDMI"
                   )


variables_used <- data.frame(variables = NULL)

StatesRankings <- list()

for (ii in seq_along(reprioritization_states)){
  
     selected_vars <- possible_vars[possible_vars %in%  names(ward_data[[ii]]) == T]
     
     # state = reprioritization_states[ii]


      filled_updata <- handle_na_neighbor_mean(data = ward_data[[ii]], 
                                               shp_data = shapefile[[ii]],
                                               col = NULL) %>% 
        # fill up the missingness in the data 
        select(WardName,WardCode, mean_EVI, mean_NDVI, mean_rainfall, distance_to_water, 
               RH_mean, temp_mean,  housing_quality, pfpr, avgRAD)
      
      

      
      normilized_data <- normalize_data(cleaned_data = filled_updata, 
                                        # why do we keeping the old datasets  
                                        variable_relationships)
      
      
      
      
      
      composite_score <- composite_score_models(normalized_data = normilized_data,
                                                # composite score for different models
                                                selected_vars = selected_vars,
                                                shp_data = shapefile[[ii]]) 
      
      
      
      plotting_data <- process_model_score(data_to_process = composite_score$final_data)
      
     
       finalranks <- plotting_data %>% 
        group_by(WardName, WardCode, flag_not_ideal) %>% 
        summarise(rank  = median(rank)) 
       
      # Note what was selected 
       variables_used <- rbind(variables_used, 
              data.frame(variables =paste("selected variables are in state",
                                   reprioritization_states[[ii]], ":",
                                   paste(selected_vars, collapse = ", "))))
       
      
      write.csv(finalranks, 
                file.path(datafile, 
                          "rankings",
                          paste0(reprioritization_states[[ii]], 
                                 "_rankings.csv")), 
                row.names = T)

}


write.csv(variables_used, 
          file.path(datafile, "rankings", "variables_per_state.csv"), 
          row.names = T)

rm(list = ls())


Drive <- gsub("\\\\", "/", Sys.getenv("USERPROFILE"))
DriveDir <- file.path(Drive, "Urban Malaria Proj Dropbox", "urban_malaria") 

shapefile <- file.path(DriveDir,"data/nigeria/nigeria_shapefiles/shapefiles/ShinyApp_shapefiles/all_reprioritization_nmep_states/")

datafile <- file.path(DriveDir,"/projects/urban_microstratification/Shiny App")



library(dplyr)
library(RColorBrewer)
library(ggiraph)
library(ggtext)


reprioritization_states <- c("Kano", "Kaduna", "Katsina", "Niger", "Taraba", "Yobe" )


handle_na_neighbor_mean <- function(data, shp_data, col = NULL) {
  
  
  if (is.null(col)) {
    cols_to_process <- names(data)[sapply(data, function(x) any(is.na(x)))]
  } else {
    cols_to_process <- col
  }
  
  # Create neighbor structure using shapefile data
  w <- spdep::poly2nb(shp_data, queen = TRUE)
  
  for (current_col in cols_to_process) {
    col_data <- data[[current_col]]
    missing_indices <- which(is.na(col_data))
    
    print(paste("Processing column:", current_col))
    print(paste("Number of NAs:", length(missing_indices)))
    
    for (index in missing_indices) {
      neighbor_indices <- w[[index]]
      neighbor_values <- col_data[neighbor_indices]
      imputed_value <- mean(neighbor_values, na.rm = TRUE)
      
      if (is.na(imputed_value)) {
        imputed_value <- mean(col_data, na.rm = TRUE)
      }
      
      col_data[index] <- imputed_value
      print(paste("Imputed value for index", index, ":", imputed_value))
    }
    
    # Only update the NAs in the original data
    data[[current_col]][missing_indices] <- col_data[missing_indices]
  }
  
  return(data)
}




# Normalization function
normalize_data <- function(cleaned_data, variable_relationships) { 
  tryCatch({
    print("Input data structure (cleaned data):")
    print(str(cleaned_data))
    print("Variable relationships:")
    print(variable_relationships) 
    
    # Identify numeric columns for normalization
    numeric_cols <- names(cleaned_data)[sapply(cleaned_data, is.numeric)]
    numeric_cols <- intersect(numeric_cols, names(variable_relationships))
    
    print("Numeric columns to be normalized:")
    print(numeric_cols)
    
    if (length(numeric_cols) == 0) {
      stop("No numeric columns found for normalization!")
    }
    
    # Apply normalization to numeric columns based on relationships
    scoring_data <- cleaned_data %>% 
      mutate(across(all_of(numeric_cols), 
                    ~{
                      col_name <- cur_column()
                      if (variable_relationships[col_name] == "inverse") {
                        inverted <- 1 / (. + 1e-10)  # Add small constant to avoid division by zero
                        ((inverted - min(inverted, na.rm = TRUE)) / 
                            (max(inverted, na.rm = TRUE) - min(inverted, na.rm = TRUE)))
                      } else {  
                        (. - min(., na.rm = TRUE)) / 
                          (max(., na.rm = TRUE) - min(., na.rm = TRUE)) 
                      }
                    },
                    .names = "normalization_{tolower(.col)}"))
    
    print("Normalized data summary:")
    print(summary(scoring_data))
    
    return(scoring_data)
    
  }, error = function(e) {
    print(paste("Error in normalize_data:", e$message))
    print(traceback()) 
    return(NULL)
  })
} 



composite_score_models <- function(normalized_data, selected_vars, shp_data) {
  
  print("Entering composite_score_models function")
  print("Normalized data structure:")
  print(str(normalized_data))
  print("Selected variables:")
  print(selected_vars)
  print("Shapefile data structure:")
  print(str(shp_data))
  
  # Get normalized column names for selected variables only
  norm_cols <- paste0("normalization_", tolower(selected_vars))
  norm_cols <- intersect(norm_cols, names(normalized_data))
  
  print("Normalized columns to be used:")
  print(norm_cols)
  
  if (length(norm_cols) < 2) {
    print("Error: At least two valid variables are required for composite score calculation.")
    return(NULL)
  }
  
  # Generate combinations
  model_combinations <- list()
  if (length(norm_cols) == 2) {
    # If only two variables are selected, create just one model
    model_combinations <- list(norm_cols)
  } else {
    # For more than two variables, generate all combinations
    for (i in 2:length(norm_cols)) {
      model_combinations <- c(model_combinations, combn(norm_cols, i, simplify = FALSE))
    }
  }
  
  if(!is.character(normalized_data$WardCode) == T){
    normalized_data$WardCode <-  as.character(normalized_data$WardCode)
  }
  
  # Calculate composite scores
  final_data <- normalized_data %>% 
    select(WardName, WardCode) %>%
    left_join(shp_data %>% select(WardName, 
                                  WardCode ,
                                  Urban), by = c("WardCode", "WardName"))
  
  for (i in seq_along(model_combinations)) {
    model_name <- paste0("model_", i)
    vars <- model_combinations[[i]]
    
    print(paste("Processing", model_name))
    print("Variables used:")
    print(vars)
    
    tryCatch({
      final_data <- final_data %>%
        mutate(!!sym(model_name) := {
          result <- rowSums(select(normalized_data, all_of(vars))) / length(vars)
          attributes(result) <- NULL # Strip attributes here
          print("Result summary:")
          print(summary(result))
          
          if (model_name == "model_4") {
            print("Detailed result for model_4:")
            print(head(result, 10))
          }
          
          # Flag if not urban and in top 5
          final_data[[paste0(model_name, "_flagged")]] <- 
            final_data$Urban == "No" & rank(result, na.last = "keep") <= 5 
          
          result
        })
    }, error = function(e) {
      print(paste("Error in composite_score_models for", model_name, ":", e$message))
      print("Data causing the error:")
      print(str(normalized_data))
      print("Variables causing the error:")
      print(vars)
    })
  }
  
  # Prepare output
  if (ncol(final_data) <= 1) {
    print("Error: No valid models could be created.")
    return(NULL)
  }
  
  list(model_formula = model_combinations, 
       final_data = final_data)
}



process_model_score <- function(data_to_process){
  # Separate Urban column
  urban_data <- data_to_process %>% select(WardName, WardCode, Urban)
  
  # Melt the data without Urban column
  melted_data <- data_to_process %>% 
    select(WardName,WardCode, starts_with("model_")) %>%  # Select model columns and WardName
    reshape2::melt(id.vars = c("WardName" = "WardName", "WardCode" = "WardCode"), variable.name = "variable", value.name = "value") 
  
  # Rejoin Urban data
  plottingdata <- melted_data %>%
    left_join(urban_data, by = c("WardName", "WardCode")) %>%
    group_by(variable) %>% 
    mutate(
      new_value = (value - min(value)) / (max(value) - min(value)),
      class = cut(new_value, seq(0, 1, 0.2), include.lowest = TRUE)
    ) %>%
    arrange(value) %>% 
    mutate(
      rank = row_number(),
      wardname_rank = paste(WardName, "(",rank,")"),
      flag_not_ideal = ifelse(Urban == "No" & rank <= 5, TRUE, FALSE)
    )
  
  print("Plotting data summary:")
  print(summary(plottingdata))
  
  plottingdata
}



plot_model_score_map <- function(shp_data, processed_csv, model_formulas, maps_per_page = 4) {
  palette_func <- brewer.pal(5, "YlOrRd")
  
  # Create facet labels with line breaks and flag
  facet_labels <- setNames(
    sapply(seq_along(model_formulas$model), function(i) {
      var_names <- strsplit(model_formulas$variables[i], " \\+ ")[[1]]
      base_label <- paste(var_names, collapse = " +<br>")
      
      # Check if the model is flagged
      if (any(processed_csv$flag_not_ideal[processed_csv$variable == model_formulas$model[i]])) {
        base_label <- paste0(base_label, "<br><span style='color:red;'>(Model Not Ideal)</span>")
      }
      
      base_label
    }),
    model_formulas$model
  )
  
  # Calculate consistent plot height based on maps per page
  plot_height <- 10 / ceiling(sqrt(maps_per_page)) 
  
  # Split the data into pages
  total_models <- nrow(model_formulas)
  pages <- ceiling(total_models / maps_per_page)
  
  plot_list <- list()
  
  for (page in 1:pages) {
    start_index <- (page - 1) * maps_per_page + 1
    end_index <- min(page * maps_per_page, total_models)
    
    current_models <- model_formulas[start_index:end_index]
    current_data <- processed_csv %>% filter(variable %in% current_models)
    
    plot <- ggplot(data = shp_data) +
      geom_sf_interactive(color = "black", fill = "white") +
      geom_sf_interactive(data = current_data, 
                          aes(geometry = geometry, fill = class, tooltip = wardname_rank)) +
      # Add a new layer for flagged wards
      geom_sf_interactive(data = current_data %>% filter(flag_not_ideal), 
                          aes(geometry = geometry), 
                          fill = NA, color = "blue", size = 1) + # Add red border
      facet_wrap(~variable, ncol = 2, labeller = labeller(variable = facet_labels)) +
      scale_fill_discrete(drop=FALSE, name="Malaria Risk Score", type = palette_func,
                          labels = c("Very Low", "Low", "Medium", "High", "Very High")) +
      labs(subtitle=paste("Page", page, "of", pages), 
           title = 'Composite Score Distribution by Model', 
           fill = "Malaria Risk Score",
           caption = "Blue outline indicates non-urban wards ranked in top 5 for de-prioritization (model may not be ideal)") +
      theme_void() +
      theme(
        strip.text = element_markdown(size = 7, face = "bold", lineheight = 1.0),
        strip.background = element_blank(), 
        legend.position = "bottom",
        legend.title = element_text(size = 6, face = "bold"),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        # Control spacing to maintain consistent plot sizes
        panel.spacing = unit(1.5, "lines"), 
        plot.caption = element_text(size = 8, hjust = 0.5)
      )
    
    # Use girafe for interactivity with fixed plot height
    plot_list[[page]] <- girafe(ggobj = plot, height_svg = plot_height) 
  }
  
  # Add a note about flagged models
  plot <- plot +
    labs(caption = "Blue outline indicates non-urban wards in top 5 (model may not be ideal)")
  
  
  return(plot_list)
}



box_plot_function <- function(plottingdata) {
  df_long <- plottingdata %>%
    select(WardName, variable, rank)
  
  ward_rankings <- df_long %>%
    group_by(WardName) %>%
    summarize(median_rank = median(rank)) %>%
    arrange(median_rank) %>%
    mutate(overall_rank = row_number())
  
  df_long <- df_long %>%
    left_join(ward_rankings, by = "WardName")
  
  df_long$WardName <- factor(df_long$WardName, levels = ward_rankings$WardName)
  
  p <- ggplot(df_long, aes(x = WardName, y = rank)) +
    geom_boxplot(fill = "#69b3a2", color = "#3c5e8b", alpha = 0.7) +
    coord_flip() +
    labs(title = "Ward Rankings Distribution", x = "", y = "Rank") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
    )
  
  plot <- ggplotly(p, height = 750) %>%
    layout(
      yaxis = list(fixedrange = FALSE),
      xaxis = list(fixedrange = TRUE)
    ) %>%
    config(scrollZoom = TRUE)
  
  return(list(plot = plot, ward_rankings = ward_rankings))
}

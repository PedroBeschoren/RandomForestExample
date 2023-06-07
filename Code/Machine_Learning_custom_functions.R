







single_physeq_to_borutaInput<-function (physeq_object, variable_to_be_classified){
  # boruta expects a transposed OTU table with variable_to_be_classified as a added variable
  # the output is a list of df ready to be used as input to boruta  
  #transpose phtseq otu table  
  otu_cells_wide_list <- #transpose the feature table...
    base::as.data.frame(t(otu_table(physeq_object)))%>%
      rownames_to_column(var = "sample")
  
  # extract sample data
  metadata_list <-
    as(sample_data(physeq_object),"data.frame")%>%
      rownames_to_column(var = "sample")
  
  #add the variable classification you want to predict with the random forest
  boruta_dataset<-
    base::merge(dplyr::select(metadata_list,sample, variable_to_be_classified),
                otu_cells_wide_list,
                by = "sample",
                all.y = TRUE)
  
  #make sure your variable to be classified is a factor, or boruta won't run
  
    boruta_dataset[,2]<-as.factor(boruta_dataset[,2]) # saves the second column, your variable_to_be_classified, as a factor
   
    output<-boruta_dataset
    
  gc()
  return(output)
}

















# write a function to split training and test data from a phyloseq object
train_and_test_spliter<-function(ps_object, variable_to_be_classified){
  
  # this function will separate train and test sets based on a phyloseq object and the variable to be predicted. it requires the function single_physeq_to_borutaInput()
  # ps_object = a phyloseq object
  # variable_to_be_classified =  a (quoted) metadata column that you want to predict
  # the output is a list of two objects: the first is the training set, the second is the test set
  
  # wrangle phyloseq data
  ps_data<-single_physeq_to_borutaInput(physeq_object = ps_object,
                                        variable_to_be_classified = variable_to_be_classified)
  
  # define training and test set. this can be ofptimized for repeated k-fold cross validation
  trainIndex<- createDataPartition(ps_data[,2], 
                                   p = .70, 
                                   list = FALSE, 
                                   times = 1)
  # set train and test sets
  data_Train <- ps_data [ trainIndex,]
  data_Test  <- ps_data [-trainIndex,]
  
  output<-list(data_Train,data_Test)
  names(output)<-c("data_Train","data_Test")
  
  return(output)
  
}





# define a function to fix borta objects and put them into formula format
fixed_boruta_formula<-function(boruta_object){
  # this fucntion takes a boruta ofbect, fixes the inconclusive tas into importnat o unimportnat, and then generates a formula
  # the input is a boruta object
  # the output is a boruta formula to be fed to caret::train
  # NOTE: boruta objects with zero imporntat features may crash!
  
  fixed_boruta<-TentativeRoughFix(boruta_object)
  boruta_imp_ASV<-getSelectedAttributes(fixed_boruta)
  print("number of importnat ASVs. Warning: if zero, formula will crash!")
  print(length(boruta_imp_ASV)%>%unlist()%>%sort())
  formula_boruta<-getConfirmedFormula(fixed_boruta)
  
  return(formula_boruta)
}


# put fixing, spliting, training and testing all in a single function
fix_split_train_test<-function (boruta_output_l, ps_object_l, variable_to_be_classified){
  # this fucntion will fix tentative features in a list of boruta objects
  # then it will split a list of phyloseq objects into training and test sets
  # then it will train the list of models
  # then it will test the lsit of models
  # then it returns a list of confusion matrixes (one for each model)
    # boruta_output_l = a list of boruta objects
    # ps_object_l = a list of phyloseq objects
    # variable_to_be_classified = the metadata varaible you are trying to predict (must be quoted, like "Stress")
  
  # NOTE ON SETTING SEED: 
  # if set.seed(3456) is kept INSIDE the function, results will be identical to another run where set.seed(3456) is kept inside the function
  # if set.seed(3456) is kept OUTSIDE the function, silencing the sed.seed inside of it, results will be identical to another run where set.seed(3456) is kept OUTSIDE the function
  # a run with set.seet(3456) INSIDE the function will be different from a run with set.seet(3456) OUTSIDE the function
  # when the function is replicated, set.seet set OUTSIDE the function will work well (different traint/test data, predictions on the replicated set; kept  consistent under the same seed for the set); time taken to calculate models will differ
  
  # fix boruta in a formula to be evaluated with caret
  boruta_formula_bac_l<-lapply(boruta_output_l, function(x) fixed_boruta_formula(x))
  
  # split train adn test dataset
  #set.seed(3456)
  train_test_l<-lapply(ps_object_l, function (x)
    train_and_test_spliter(ps_object = x, 
                           variable_to_be_classified = variable_to_be_classified))
  
  
  
  # train model
  #set.seed(3456)
  boruta_feature_rf_repeatedcv<-mapply(function (x,z) {
    
    train.control <- caret::trainControl(method = "repeatedcv", # set trainig/data split controls for the train function
                                         number = 5,
                                         repeats = 50,
                                         allowParallel = TRUE)
    
    model_borutized <- caret::train(form = z, # bruta formula
                                    data = x[[1]], # training data ; first element of train_and_test_spliter()
                                    method = "rf", #execute training based on RF
                                    trControl = train.control, # defined in trainControl() above
                                    ntree=2500)
    
    
    
    return(model_borutized)
  },
  x = train_test_l,
  z = boruta_formula_bac_l,
  SIMPLIFY = FALSE)
  
  
  
  #test model
 # set.seed(3456)
  confusion_matrix_output<-mapply(function(x,y){
    prediction<-stats::predict(object = x, newdata = y[[2]]) 
    confusion_output<-confusionMatrix(data = prediction, reference = y[[2]][,2])
    return(confusion_output)
  },
  x = boruta_feature_rf_repeatedcv,
  y = train_test_l,
  SIMPLIFY = FALSE)
  
  output<-list("trained_model_rf_repeatedcv" = boruta_feature_rf_repeatedcv,
               "confusion_matrix_output" = confusion_matrix_output)
  
  return(output)
  
  
}




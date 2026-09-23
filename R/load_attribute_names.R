#' Load an attributes file (.csv or .dbf) and get attribute names
#'
#' @param attributes_path Path to attributes file. Preferable to use the RangeMap_Attribute.csv so full field names are preserved, but the .tif.vat.dbf associated with one year's raster also works.
#' @return Names of all numeric attributes that can be generated
#' @export
#'
load_attribute_names<- function(attributes_path){
  
  # Check if the attributes path is csv or dbf, then load it
  if(endsWith(attributes_path, ".dbf")){
    atts<- foreign::read.dbf(attributes_path)
    message("Attributes are from dbf file. Recommend using RangeMap Attributes csv for full attribute names")
  } else if(endsWith(attributes_path, ".csv")){
    atts<- utils::read.csv(attributes_path, check.names = F)
  } else{
    message("Expected .csv or .dbf file, something else provided")
  }
  
  # Get attribute names, and remove extraneous data
  attribute_names<- names(atts)
  attribute_names<- attribute_names[-which(attribute_names %in% c("Value", "RM_ID", "Count", "PrimaryKey", "DataSource"))]
  #
  return(attribute_names)
}



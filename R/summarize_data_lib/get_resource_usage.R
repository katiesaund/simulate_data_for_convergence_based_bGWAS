# ------------------------------------------------------------------------------
# GetResourceUsage
# 2020/01/29
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Take a job id and return the Job Wall-clock time in hours and the Memory 
# Utilized in GB as a numeric vector
# ------------------------------------------------------------------------------

GetResourceUsage = function( jobId )
{
  # Run the system command to get the seff result for this job ID
  seff = system( paste("seff", jobId ), intern = TRUE)
  
  # Find the position of the time and the memory utilized in the 
  # seff result
  wallClockStr = "Job Wall-clock time: "
  wallTime = gsub(wallClockStr, '', seff[ grep(wallClockStr, seff) ])
  memoryUseStr   = "Memory Utilized: "
  memoryUsage = gsub(memoryUseStr, '', seff[ grep(memoryUseStr, seff) ])
  
  # If memory usage from seff if returned in MB, convert to GB
  if ( grepl("MB", memoryUsage) )
  {
    memoryUsage = gsub("MB", '', memoryUsage)
    memoryUsage = as.numeric( trimws(memoryUsage) ) / 1000
  } else if ( grepl("GB", memoryUsage) ) {
    memoryUsage = gsub("GB", '', memoryUsage)
    memoryUsage = as.numeric( trimws(memoryUsage) )
  } else {
    memoryUsage = gsub("KB", '', memoryUsage)
    memoryUsage = as.numeric( trimws(memoryUsage) ) * 1000
  }
  
  return( c(ConvertTimeToHours(wallTime), memoryUsage) )
}

# Take the seff time result in hours:minutes:seconds and convert the 
# results to hours
ConvertTimeToHours = function( timeStr )
{
  times   = as.numeric( unlist(strsplit(timeStr, ':') ) )
  expVals = c(0, 1, 2)
  hours   = sum( sapply( seq(3), function(i) times[i] / 60**expVals[i] ) )
  return(  round(hours, 3) )
}

# ------------------------------------------------------------------------------

#'Generate size distribution of sunflower pyrite based on different growth mechanisms
#'
#' @param Number_Sunflowers An integer indicanting the number of sunflowers to grow
#' @param Initialization If TRUE the growth starts from a framboids data_frame.
#'    If FALSE the growth continues from a sunflower data_frame
#' @param framboids When Initialization=TRUE contains the framboid data_frame from which sunflower grow
#' @param sunflowers When Initialization=FALSE contains a sunflower data_frame previously grown
#' @param Simulation The growth mechanism to simulate.
#'    1,3 surface-controlled,
#'    2-4 supply-controlled,
#'    1,2 size dependent,
#'    3,4 size independent,
#'    5 supply-controlled, increasing volume, size dependent
#' @param Random_Limit A number from 0 to 1 representing the maximum value of random numbers
#'    in algorithm
#' @param Iterations An integer value that controls the maximum number of iterations of the algorithm
#' @param MaxInfillingVolume A value from 0 to 1-packed_factor indicating the maximim infilled volume
#'  of the framboidal core (Simulation = 5)
#' @param MaxMeanDiameter A number (in micrometers) controling the maximum mean of the size
#'    distribution of sunflowers
#' @return A sunflower data_frame
#' @seealso \code{\link{Framboids_growth}}
#' #' @references
#' ## Merinero, R.; Cardenes, V. (in press)
#' ## Theoretical growth of framboidal and sunflower pyrite using the R-package frambgrowth.
#' ## Mineralogy and Petrology. doi:10.1007/s00710-017-0535-x
#' @examples
#' ## Size dependent and surface-controlled growth of 100 framboids adding nanocrystals
#' ## until the mean of the diameters of framboids was 10 micrometers followed by size dependent
#' ## and supply-controlled growth until the mean of the diameter was 15 micrometers followed by
#' ## supply-controlled, increasing volume and size dependent growth of sunflowers
#' ## until the mean of the diameter was 20 micrometers
#' library(frambgrowth)
#' Framboids<-Framboids_growth(100,Initialization=TRUE, Simulation=1,
#' Iterations=1000, MaxMeanDiameter=10)
#'Framboids2<-Framboids_growth(100,Initialization=FALSE, framboids=Framboids, Simulation=5,
#' Iterations=1000, MaxMeanDiameter=15)
#' Sunflowers<-Sunflowers_growth(100,Initialization=TRUE, framboids=Framboids2, Simulation=5,
#' Iterations=1000, MaxMeanDiameter=20)


Sunflowers_growth <- function(Number_Sunflowers=1000, Initialization=FALSE, framboids, sunflowers,
                             Simulation=1, Random_Limit=1, Iterations=3, MaxInfillingVolume=0.1,
                             MaxMeanDiameter=20)
{
  ## Initial testing
  if (Number_Sunflowers<1) {print('Number of sunflowers lower than 0'); return()}
  if ((!(Initialization))&&(is.na(sunflowers))) {print('Please specify the variable sunflowers');return()}

  ## Auxiliars and work variables
  N<- Number_Sunflowers
  mu <- Random_Limit
  MaxVr <- MaxInfillingVolume

  ## sunflowers data.frame contains the growing sunflowers and it is the return of the function
  ## Initialization of the data.frame if it is indicated in parameters

  if (Initialization)
  {
    sunflowers<-data.frame(c = 1:N, diameter = 0, volume = 0, external_volume = 0, infilled_volume = 0,
                           framboidal_diameter = 0, framboidal_initial_diameter = 0, framboidal_nanocrystals = 0,
                           framboidal_nanocrystals_diameter =0, framboidal_packed_factor = 0, framboidal_volume =0)
    for(j in 1:N)
    {
    sunflowers$framboidal_diameter <- framboids$diameter
    sunflowers$diameter <- sunflowers$framboidal_diameter
    sunflowers$framboidal_initial_diameter <- framboids$initial_diameter
    sunflowers$framboidal_nanocrystals <- framboids$nanocrystals
    sunflowers$framboidal_nanocrystals_diameter <- framboids$nanocrystals_diameter
    sunflowers$framboidal_packing_factor <- framboids$nanocrystals_factor
    sunflowers$framboidal_volume <- framboids$volume
    sunflowers$volume <- sunflowers$framboidal_volume
    }
  }

  #### Random initialization
  runif(1)

  #### i = growth cycles
  i<-0
  ####

  ## Average diameter of the framboids is one of the controls for the loop

  AverageDiameter <- mean(sunflowers$diameter)

  ## Loop of growth cycles controled by the number of iterations passed to the function or the MaxMeanDiameter expected

  while ((i < Iterations) && (AverageDiameter < MaxMeanDiameter))
  {

    ## EpsilonC random variable of the growth cycle
    EpsilonC <- runif(1)*mu

    ## j = sunflowers

    for(j in 1:N)
    {

      ## EpsilonF random variable for each framboid j in each growth cycle i

      EpsilonF <- runif(1)*mu

      ## Initialization of empty sunflowers

      if (sunflowers$diameter[j]==0)
      {
        sunflowers$diameter[j]<-round(1+EpsilonF)
      }

      else

      {

        if (Simulation == 1)
        {

          ### Simulation = 1 surface-controlled, increasing diameter and size dependent

          sunflowers$diameter[j] <- sunflowers$diameter[j]*(1+EpsilonF)

        }

        if (Simulation == 2)
        {

          ### Simulation = 2 supply-controlled, increasing diameter and size dependent

          sunflowers$diameter[j] <- sunflowers$diameter[j]*(1+EpsilonF*EpsilonC)

        }

        if (Simulation == 3)
        {

          ### Simulation = 3 surface-controlled, increasing diameter and size independent

          sunflowers$diameter[j] <- sunflowers$diameter[j] + EpsilonF

        }

        if (Simulation == 4)
        {

          ### Simulation = 4 supply-controlled, increasing diameter and size independent

          sunflowers$diameter[j] <- sunflowers$diameter[j] + EpsilonF*EpsilonC

        }

        if (Simulation == 5)
        {

          ### Simulation = 5 supply-controlled, increasing volume, size dependent
          ### with previous infilling of the framboidal core (if indicated)

          landavolume<-sunflowers$volume[j]*EpsilonF*EpsilonC
          F <- sunflowers$framboidal_packed_factor[j]
          D0 <- sunflowers$framboidal_initial_diameter[j]
          W<- ((4*pi*F)/3)*((D0/2)^3)
          volumethframboid<-((4*pi)/3)*((sunflowers$framboidal_diameter[j]/2)^3) - W
          if (volumethframboid*MaxVr > landavolume + sunflowers$infilled_volume[j])
          {
            ### Volume added to the infilling of the framboidal core

            sunflowers$infilled_volume[j]<-sunflowers$infilled_volume[j]+landavolume
          }
          else

            ### Volume added to the external regrowth

          {
            sunflowers$external_volume[j]<-sunflowers$external_volume[j]+landavolume
            + sunflowers$infilled_volume[j] - volumethframboid*MaxVr
            sunflowers$infilled_volume[j]<-volumethframboid*MaxVr
          }

          ### The diameter is recalculated as function of the volume

          sunflowers$diameter[j] <- 2*((3*(volumethframboid+sunflowers$external_volume[j]))/(4*pi))^(1/3)

        }
      }

      ## External volume and volume of the sunflower

      if (Simulation != 5)
      {
        sunflowers$external_volume[j] <- ((4*pi)/3)*((sunflowers$diameter[j]/2)^3) - ((4*pi)/3)*((sunflowers$framboidal_diameter[j]/2)^3)
        sunflowers$volume[j] <- sunflowers$external_volume[j] + sunflowers$framboidal_volume[j] + sunflowers$infilled_volume[j]
      }


    }

    i = i + 1
    AverageDiameter <- mean(sunflowers$diameter)

  }
return(sunflowers)
}


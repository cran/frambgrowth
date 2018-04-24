#' Generate size distribution of framboidal pyrite based on different growth mechanisms
#'
#' @param Number_Framboids An integer indicanting the number of framboids to grow
#' @param Initialization If TRUE the growth starts from empty framboids. If FALSE the growth continues
#'    from framboids
#' @param framboids When Initialization=FALSE contains a framboids data_frame previously grown
#' @param Simulation The growth mechanism to simulate.
#'    1-4 surface-controlled,
#'    5-6 supply-controlled,
#'    1,3,5,6 size dependent,
#'    2,4,7,8 size independent,
#'    1,2,5,7 adding nanocrystals,
#'    3,4,6,8 increasing diameter
#' @param Initial_Diameter A number (in micrometers) representing the growth over preexistent
#'    spherical objects
#' @param Nanocrystals_Diameter A number (in micrometers) with the diameter of nanocrystals
#'    forming framboids
#' @param Packing_Factor A number from 0 to 1 representing the packing factor of nanocrystals
#' @param Random_Limit A number from 0 to 1 representing the maximum value of random numbers
#'    in algorithm
#' @param Iterations An integer value that controls the maximum number of iterations of the algorithm
#' @param MaxMeanDiameter A number (in micrometers) controlling the maximum mean of the size
#'    distribution of framboids
#' @return A framboids data_frame
#' @seealso \code{\link{Sunflowers_growth}}
#' @references
#' ## Merinero, R.; Cardenes, V. (in press)
#' ## Theoretical growth of framboidal and sunflower pyrite using the R-package frambgrowth.
#' ## Mineralogy and Petrology. doi:10.1007/s00710-017-0535-x
#' @examples
#' ## Size dependent and surface-controlled growth of 1000 framboids adding nanocrystals
#' ## until the mean of the diameters of framboids was 15 micrometers.
#' library(frambgrowth)
#' Framboids<-Framboids_growth(100,Initialization=TRUE, Simulation=1,
#' Iterations=1000, MaxMeanDiameter=15)
#' ## Size dependent and surface-controlled growth of 1000 framboids adding nanocrystals
#' ## until the mean of the diameters of framboids was 10 micrometers,
#' ## followed by size dependent and supply-controlled growth until the mean was 15 micrometers.
#' Framboids<-Framboids_growth(100,Initialization=TRUE, Simulation=1,
#' Iterations=1000, MaxMeanDiameter=10)
#' Framboids<-Framboids_growth(100,Initialization=FALSE, Simulation=5, framboids=Framboids,
#' Iterations=1000, MaxMeanDiameter=15)



Framboids_growth <- function(Number_Framboids=1000,Initialization=TRUE,framboids,Simulation=1,
                             Initial_Diameter=0, Nanocrystals_Diameter=0.1,Packing_Factor=0.74078,
                             Random_Limit=1, Iterations=3, MaxMeanDiameter=20)




{

  ## Initial testing
  if (Number_Framboids<1) {print('Number of framboids lower than 0'); return()}
  if ((!(Initialization))&&(is.na(framboids))) {print('Please specify the variable framboids');return()}

  ## Auxiliars and work variables

  N<- Number_Framboids
  mu <- Random_Limit

  ## framboids data.frame contains the growing framboids and it is the return of the function
  ## Initialization of the data.frame if it is indicated in parameters

  if (Initialization)
  {
    D0<- Initial_Diameter
    d<- Nanocrystals_Diameter
    F<- Packing_Factor

    framboids<-data.frame(c = 1:N, diameter = D0, initial_diameter = D0, nanocrystals_diameter = d,
                          nanocrystals_factor = F, volume = 0, nanocrystals = 0)

    ## Initial diameter and nanocrystals diameter can be sligthly different between framboids

    for(j in 1:N)
    {

      if (D0 > 0)
      {
        framboids$initial_diameter[j] <- D0 + runif(1)*0.01
        framboids$diameter[j] <- framboids$initial_diameter[j]
      }

      framboids$nanocrystals_diameter[j] <- (d-0.01) + runif(1)*0.02
      framboids$nanocrystals_factor[j] <- (F-0.01) + runif(1)*0.02

    }

  }

    #### Random initialization
    runif(1)

    #### i = growth cycles
    i<-0
    ####

    ## Average diameter of the framboids is one of the controls for the loop

    AverageDiameter <- mean(framboids$diameter)

    ## Loop of growth cycles controled by the number of iterations passed to the function or the MaxMeanDiameter expected

    while ((i < Iterations) && (AverageDiameter < MaxMeanDiameter))
    {

      ## EpsilonC random variable of the growth cycle
      EpsilonC <- runif(1)*mu

      ## j = framboids

      for(j in 1:N)
      {

        ## EpsilonF random variable for each framboid j in each growth cycle i

        EpsilonF <- runif(1)*mu

        ## Work variables
        F <- framboids$nanocrystals_factor[j]
        d <- framboids$nanocrystals_diameter[j]
        D0 <- framboids$initial_diameter[j]
        W<- ((4*pi*F)/3)*((D0/2)^3)
        P<- round(F*((D0/d)^3))

        ## Initialization of empty framboids
        if (framboids$nanocrystals[j]==0)
        {
          framboids$nanocrystals[j]<-round(1+EpsilonF)
        }

        else

        {

          if (Simulation == 1)
          {
            ### Simulation = 1 surface-controlled, adding nanocrystals and size dependent

            framboids$nanocrystals[j]<-framboids$nanocrystals[j]*(1+EpsilonF)

          }

          if (Simulation == 2)

          {
            ### Simulation = 2 surface-controlled, adding nanocrystals and size independent

            framboids$nanocrystals[j]<-framboids$nanocrystals[j] + EpsilonF

          }

          if (Simulation == 3)
          {
            ### Simulation = 3 surface-controlled, incresing diameter and size dependent

            framboids$diameter[j]<-framboids$diameter[j]*(1+EpsilonF)
            framboids$nanocrystals[j]<-round(F*((framboids$diameter[j]/d)^3)-P)

          }


          if (Simulation == 4)

          {
            ### Simulation = 4 surface controlled, increasing diameter and size independent

            framboids$diameter[j]<-framboids$diameter[j] + EpsilonF
            framboids$nanocrystals[j]<-round(F*((framboids$diameter[j]/d)^3)-P)

          }

          if (Simulation == 5)

          {
            ### Simulation = 5 supply controlled, adding nanocrystals and size dependent

            framboids$nanocrystals[j]<-framboids$nanocrystals[j]*(1+EpsilonF*EpsilonC)

          }

          if (Simulation == 6)

          {
            ### Simulation = 6 supply-controled, increasing diameter and size dependent

            framboids$diameter[j]<-framboids$diameter[j]*(1+EpsilonF*EpsilonC)
            framboids$nanocrystals[j]<-round(F*((framboids$diameter[j]/d)^3)-P)

          }

          if (Simulation == 7)

          {
            ### Simulation = 7 supply controlled, adding nanocrystals and not size dependent

            framboids$nanocrystals[j]<-framboids$nanocrystals[j] + EpsilonF*EpsilonC

          }

          if (Simulation == 8)

          {
            ### Simulation = 8 supply-controled, increasing diameter and not size dependent

            framboids$diameter[j]<-framboids$diameter[j] + EpsilonF*EpsilonC
            framboids$nanocrystals[j]<-round(F*((framboids$diameter[j]/d)^3)-P)

          }
        }

        ## Diameter and volume calculation

        framboids$diameter[j]<- d*(((framboids$nanocrystals[j]+P)/F)^(1/3))
        framboids$volume[j] <- ((F*4*pi)/3)*(( framboids$diameter[j]/2)^3)-W

      }
      ### print(c(mean(framboids$diameter),sd(framboids$diameter),max(framboids$diameter)))

      i <- i + 1
      AverageDiameter <- mean(framboids$diameter)


    }
    return(framboids)
  }

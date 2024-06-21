# TODOs:
# - implement treatment of B as T1
# - normalize measurements, show panels per subject
# - note: date of birth is available in the CovidData, this can't be published!

# ------------------------------------------------------------------------

ImmunData <- R6::R6Class(
  "ImmunData",
  private = list(
    rawData = NULL
  ),
  public = list(
    initialize = function( fileName ) {
      rawData <- readxl::read_excel( fileName ) |>
        dplyr::rename( SubjectId = study_ID, SampleId = Sample_ID ) |>
        tidyr::drop_na( SampleId ) |>
        dplyr::mutate( RowId = sprintf( "ROW_%05d", dplyr::row_number() ) )
      
      isDup <- duplicated( rawData$SampleId )
      if( any( isDup ) ) {
        warning( "Removing rows for duplicated SampleId: ", rawData$SampleId[ isDup ] )
        rawData <- rawData[ !isDup, ]
      }
      stopifnot( all( !duplicated( rawData$SampleId ) ) )
      private$rawData <- rawData
    },
    getData = function() {
      data <- private$rawData
      data
    },
    getMeasuresData = function() {
      data <- self$getData()
      md <- data |>
        tidyr::drop_na( SubjectId ) |>
        dplyr::select( -sample_id, -Data_measured, -RowId ) |>
        tidyr::pivot_longer(
          cols = c( -SubjectId, -SampleId ),
          names_to = "Marker",
          values_to = "Value"
        ) |>
        dplyr::arrange( Marker, SubjectId, SampleId ) |>
        tidyr::drop_na( Value )
      keepMarkers <- md |>
        dplyr::group_by( Marker ) |>
        dplyr::summarize( keep = max( Value ) != min( Value ), .groups = "drop" ) |>
        dplyr::filter( keep ) |>
        dplyr::pull( Marker )
      md <- md |>
        dplyr::filter( Marker %in% keepMarkers ) |>
        dplyr::select( Marker, SubjectId, SampleId, Value, dplyr::everything() )
      md
    },
    getMarkers = function() {
      self$getMeasuresData() |>
        dplyr::mutate( dup = duplicated( Marker ) ) |>
        dplyr::filter( !dup ) |>
        dplyr::pull( Marker )
    }
  )
)

if( FALSE ) {
  immData <- ImmunData$new( "Dataset/IMMUN_COVID20240130_V2.xlsx" )
  d <- immData$getData(); d
  md <- immData$getMeasuresData(); md
}

# ------------------------------------------------------------------------

CovidData <- R6::R6Class(
  "CovidData",
  private = list(
    rawData = NULL,
    vaxValidDaysRange = NULL,
    reinfDaysRange = NULL
  ),
  public = list(
    initialize = function( 
      fileName, 
      vaxValidDaysRange = c( 14, 270 ), reinfDaysRange = c( 14, 270 )
    ) {
      private$vaxValidDaysRange <- vaxValidDaysRange
      private$reinfDaysRange <- reinfDaysRange

      rawData <- readxl::read_excel( fileName ) |>
        dplyr::mutate( RowId = sprintf( "ROW_%05d", dplyr::row_number() ) ) |>
        dplyr::rename( SampleId = Sample_ID, SubjectId = ZMW_ID )

      rawData <- rawData |> 
        dplyr::mutate(
          CollectionDate = lubridate::dmy(CollectionDate),
  #        COVID1 = lubridate::dmy(COVID1),
  #        COVID2 = lubridate::dmy(COVID2),
  #        COVID3 = lubridate::dmy(COVID3),
  #        COVID4 = lubridate::dmy(COVID4),
          Vax1 = lubridate::dmy(Vax1),
          Vax2 = lubridate::dmy(Vax2),
          Vax3 = lubridate::dmy(Vax3),
          Vax4 = lubridate::dmy(Vax4),
          Vax5 = lubridate::dmy(Vax5),
          Vax6 = lubridate::dmy(Vax6)
          #,        Project = str_replace_all(Project, c(" " = "", "-" = "_"))
        )

      private$rawData <- rawData
    },
    getData = function() {
      data <- private$rawData
      data
    },
    getSubjectData = function() {
      data <- self$getData()

      subData <- data |>
        dplyr::distinct( SubjectId, Date_of_birth, Gender )
      stopifnot( !any( duplicated( subData$SubjectId ) ) )
      subData
    },
    getVaccinationsData = function() {
      data <- self$getData()

      vd <- data |>
        dplyr::select(
          SubjectId,
          Vax1, Vax2, Vax3, Vax4, Vax5, Vax6
        ) |>
        tidyr::pivot_longer(
          cols = c( -SubjectId ),
          names_to = "Event",
          values_to = "Date"
        ) |>
        tidyr::drop_na( Date ) |>
        dplyr::distinct( SubjectId, Date ) |>
        dplyr::arrange( SubjectId, Date )
      vd
    },
    getCollectionData = function() {
      data <- self$getData()
      vd <- self$getVaccinationsData()

      cd <- data |>
        dplyr::select( SampleId, SubjectId, Date = CollectionDate, CaseVSprevious, Case_type ) |>
        tidyr::drop_na( Date )
      infD <- cd |>
        dplyr::filter( CaseVSprevious %in% c( "CASE" ) & Case_type %in% c( "asympt", "sympt" ) )

      # ----- calculating vaccination breakthrough status -----
      vaxValidDaysRange <- private$vaxValidDaysRange
      d <- infD |>
        dplyr::inner_join( 
          vd |> dplyr::rename( VaxDate = Date ),
          dplyr::join_by( SubjectId == SubjectId, dplyr::closest( Date > VaxDate ) )
        ) |>
        dplyr::mutate( 
          isVaxBT = 
            ( Date - VaxDate ) >= lubridate::days( 14 + min( vaxValidDaysRange ) ) & 
            ( Date - VaxDate ) <= lubridate::days( 14 + max( vaxValidDaysRange ) )
        ) |>
        dplyr::select( SampleId, SubjectId, VaxDate, isVaxBT )
      cd <- cd |>
        dplyr::left_join( d, by = c( "SampleId", "SubjectId" ) )

      # ----- calculating reinfection status -----
      reinfDaysRange <- private$reinfDaysRange
      d <- infD |>
        dplyr::inner_join( 
          infD |> dplyr::select( SubjectId, InfDate = Date ),
          dplyr::join_by( SubjectId == SubjectId, dplyr::closest( Date > InfDate ) )
        ) |>
        dplyr::mutate( 
          isReinf = 
            ( Date - InfDate ) >= lubridate::days( min( reinfDaysRange ) ) & 
            ( Date - InfDate ) <= lubridate::days( max( reinfDaysRange ) )
        ) |>
        dplyr::select( SampleId, SubjectId, InfDate, isReinf )
      cd <- cd |>
        dplyr::left_join( d, by = c( "SampleId", "SubjectId" ) )

      cd
    },
    getMeasuresData = function() {
      data <- self$getData()
      
      md <- data |>
        dplyr::select(
          SampleId, SubjectId, Date = CollectionDate,
          BA1, BA2, BA5, NCP, S1, QFNS_Ag1, QFNS_Ag2
        ) |>
        tidyr::pivot_longer(
          cols = c( -SampleId, -SubjectId, -Date ),
          names_to = "Marker",
          values_to = "Value"
        ) |>
        tidyr::drop_na()

      md <- md |>
        dplyr::arrange( Marker, SubjectId, Date ) |>
        dplyr::select( Marker, SubjectId, Date, Value, SampleId, dplyr::everything() )

      md
    },
    getMarkers = function() {
      self$getMeasuresData() |>
        dplyr::mutate( dup = duplicated( Marker ) ) |>
        dplyr::filter( !dup ) |>
        dplyr::pull( Marker )
    }
  )
)

if( FALSE ) {
  covData <- CovidData$new( "Dataset/BREAK COVID overview v3_combined_data__20240102_LUMC_Mo_SzMK.xlsx" )
  data <- covData$getData(); data
  subData <- covData$getSubjectData(); subData
  md <- covData$getMeasuresData(); md
  vd <- covData$getVaccinationsData(); vd
  cd <- covData$getCollectionData(); cd
}

# ------------------------------------------------------------------------

DataMerger <- R6::R6Class(
  "DataMerger",
  private = list(
    transDaysRange = NULL,
    preferredDaysDiff = NULL,
    immunData = NULL,
    covidData = NULL,
    projectSources = NULL
  ),
  public = list(
    initialize = function( 
      immunData, covidData, 
      transDaysRange = c( 14, 9*31 ), preferredDaysDiff = 180, projectSources = NULL
    ) {
      private$transDaysRange <- transDaysRange
      private$preferredDaysDiff <- preferredDaysDiff
      private$immunData <- immunData
      private$covidData <- covidData
      private$projectSources <- projectSources
    },
    getVaccinationsData = function() {
      vd <- private$covidData$getVaccinationsData()
      vd
    },
    getMeasures = function( eventType = "all" ) {
      transDaysRange <- private$transDaysRange
      preferredDaysDiff <- private$preferredDaysDiff
      covData <- private$covidData
      immData <- private$immunData
      projectSources <- private$projectSources

      # ----- merge meassurements -----
      md <- dplyr::bind_rows(
        immData$getMeasuresData() |> dplyr::select( Marker, SampleId, Value ) |> dplyr::mutate( Source = "imm" ),
        covData$getMeasuresData() |> dplyr::select( Marker, SampleId, Value ) |> dplyr::mutate( Source = "cov" )
      )

      if( !is.null( projectSources ) ) {
        d <- covData$getData() |> dplyr::select( SampleId, ProjectSource ) |> dplyr::filter( ProjectSource %in% projectSources )
        md <- md |>
          dplyr::semi_join( d, by = "SampleId" )
      }

      # ----- add dates of collection -----
      md <- md |>
        dplyr::left_join( covData$getCollectionData(), by = "SampleId" )

      # ----- calculate age at date of collection -----
      subData <- covData$getSubjectData()
      md <- md |>
        dplyr::left_join( subData, by = "SubjectId" ) |>
        dplyr::mutate( AgeYears = as.numeric( Date - lubridate::dmy( Date_of_birth ), "days" ) / 365.25 ) |>
        dplyr::select( -Date_of_birth )

      # ----- number measurements individually for each marker/subject -----
      md <- md |>
        dplyr::group_by( Marker, SubjectId ) |>
        dplyr::mutate( DatePoint = rank( Date, ties.method = "random" ) ) |>
        dplyr::ungroup() |>
        dplyr::mutate( DatePoint = ordered( DatePoint ) )

      if( eventType == "all" ) {
        mInfo <- tibble::tribble(
          ~CaseVSprevious, ~Case_type, ~Event, ~EventCode, ~EventColor, ~EventState,
          "baseline", NA, "Baseline", "B", "darkolivegreen", "control",
          "CONTROL", NA, "Control", "C", "darkgreen", "control",
          "CASE", "asympt", "Asymptomatic", "A", "orange", "case",
          "CASE", "sympt", "Symptomatic", "S", "red", "case",
          "CASE", "dubious", "Dubious", "D", "gray", NA
        )
        md <- md |>
          dplyr::left_join( mInfo |> dplyr::select( -EventColor ), by = c( "CaseVSprevious", "Case_type" ) ) |>
          dplyr::mutate( Event = factor( Event, levels = mInfo$Event ) ) |>
          dplyr::mutate( EventState = factor( EventState, levels = c( "control", "case" ) ) )
        mInfo <- mInfo |> dplyr::select( Event, EventCode, EventColor )
      } else {
        stop( "Unknown eventType: ", eventType )
      }

      # ----- preparing vax/reinf codes -----
      if( eventType == "all" ) {
        vrInfo <- tibble::tribble(
          ~VaxReinfCode, ~VaxReinfColor,
          "VBT/ReInf", "purple",
          "VBT", "blue",
          "ReInf", "red",
          ".", "darkgrey"
        )
        vrMap <- tibble::tribble(
          ~isVaxBT, ~isReinf, ~VaxReinfCode,
          TRUE, TRUE, "VBT/ReInf",
          TRUE, FALSE, "VBT",
          TRUE, NA, "VBT",
          FALSE, TRUE, "ReInf",
          NA, TRUE, "ReInf",
          FALSE, FALSE, ".",
          FALSE, NA, ".",
          NA, FALSE, ".",
          NA, NA, "."
        )
        md <- md |>
          dplyr::left_join( vrMap, by = c( "isVaxBT", "isReinf" ) ) |>
          dplyr::mutate( VaxReinfCode = factor( VaxReinfCode, levels = vrInfo$VaxReinfCode ) )
      } else {
        stop( "Unknown eventType: ", eventType )
      }

      # ----- preparing transitions -----
      td <- md |>
        dplyr::select( Marker, SubjectId, SampleId, Date, Value, Event, EventCode )
      prevTd <- td |> 
        dplyr::rename( prevDate = Date, prevValue = Value, prevEvent = Event, prevEventCode = EventCode, prevSampleId = SampleId )
      td <- td |>
#        dplyr::inner_join( prevTd, dplyr::join_by( SubjectId, Marker, dplyr::closest( Date > prevDate ) ) )
        dplyr::inner_join( prevTd, dplyr::join_by( SubjectId, Marker, Date > prevDate ) ) |>
        dplyr::arrange( SubjectId, Marker, prevDate, Date )
      td <- td |>
        dplyr::filter( 
          ( Date - prevDate ) >= lubridate::days( min( transDaysRange ) ) & 
          ( Date - prevDate ) <= lubridate::days( max( transDaysRange ) )
        )

      if( eventType == "all" ) {
        ttName2rank <- c( "BA" = 1L, "BS" = 1L, "CA" = 1L, "CS" = 1L, "CC" = 2L, "BC" = 2L )
        tInfo <- dplyr::bind_rows(
          tidyr::expand_grid( pe = c( "B", "C" ), e = c( "B", "C" ), c = "darkgreen", n = "control.control" ),
          tidyr::expand_grid( pe = c( "B", "C" ), e = c( "S", "A" ), c = "red", n = "control.case" ),
          tidyr::expand_grid( pe = c( "S", "A" ), e = c( "B", "C" ), c = "lightblue", n = "case.control" ),
          tidyr::expand_grid( pe = c( "S", "A" ), e = c( "S", "A" ), c = "bisque3", n = "case.case" ),
          tidyr::expand_grid( pe = c( "B", "C" ), e = c( "D" ), c = "gray", n = NA ),
          tidyr::expand_grid( pe = c( "S", "A" ), e = c( "D" ), c = "gray", n = NA ),
          tidyr::expand_grid( pe = c( "D" ), e = c( "B", "C", "A", "S", "D" ), c = "gray", n = NA ),
        ) |>
          dplyr::mutate( TransitionRank = ttName2rank[ paste0( pe, e ) ] ) |>
          dplyr::rename( prevEventCode = pe, EventCode = e, TransitionName = n, TransitionColor = c )

        td <- td |>
          dplyr::left_join( 
            tInfo |> dplyr::select( -TransitionColor ),
            by = c( "prevEventCode", "EventCode" )
          ) |>
          tidyr::drop_na( TransitionName ) |>
          tidyr::drop_na( TransitionRank ) |>
          dplyr::mutate( DateDiff = abs( preferredDaysDiff - as.integer( Date - prevDate ) ) ) |>
          dplyr::group_by( Marker, SubjectId ) |>
          dplyr::filter( TransitionRank == min( TransitionRank ) ) |>
          dplyr::slice_min( DateDiff, n = 1L, with_ties = FALSE ) |>
          dplyr::ungroup()
      }

      # ----- result -----
      list( 
        measurementData = md, measurementInfo = mInfo, vaxReinfInfo = vrInfo,
        transitionData = td, 
        transitionInfo = tInfo |> tidyr::drop_na( TransitionName )
      )
    }
  )
)

if( FALSE ) {
  covData <- CovidData$new( "Dataset/BREAK COVID overview v3_combined_data__20240102_LUMC_Mo_SzMK.xlsx" )
  immData <- ImmunData$new( "Dataset/IMMUN_COVID20240130_V2.xlsx" )
  dm <- DataMerger$new( immData, covData )
  mm <- dm$getMeasures(); mm
  vd <- dm$getVaccinationsData(); vd
}

# ------------------------------------------------------------------------

DataProcessor <- R6::R6Class(
  "DataProcessor",
  private = list(
    dataMerger = NULL
  ),
  public = list(
    initialize = function( dataMerger ) {
      private$dataMerger <- dataMerger
    },
    summary = function() {
      dm <- private$dataMerger
      mm <- dm$getMeasures() 
      mmDatesData <- mm$data |>
        dplyr::distinct( Date, SubjectId, EventCode, isVaxBT, isReinf )

      l <- list(
        vaxBTreinfInfo = mmDatesData |> dplyr::count( isVaxBT, isReinf )
      )
      l
    },

    # ----- data providers -----
    getVaccinationsData = function() {
      dm <- private$dataMerger
      dm$getVaccinationsData()
    },
    getMeasures = function() {
      dm <- private$dataMerger
      dm$getMeasures() 
    },
    getSubjectIds = function() {
      mm <- self$getMeasures()
      sIds <- mm$measurementData$SubjectId
      sIds <- sIds[ !duplicated( sIds ) ]
      sIds
    },
    getMarker2measurementData = function( markers = NULL ) {
      mm <- self$getMeasures() 
      m2md <- split( mm$measurementData, mm$measurementData$Marker )

      if( !is.null( markers ) ) {
        m2md <- m2md[ intersect( markers, names( m2md ) ) ]
      }

      m2md
    },
    getMarker2transitionData = function( twoLevels = TRUE, markers = NULL ) {
      mm <- self$getMeasures() 
      m2td <- split( mm$transitionData, mm$transitionData$Marker )

      if( twoLevels ) {
        m2td <- lapply( m2td, function( td ) {
          td |>
            dplyr::mutate( TransitionName = factor( TransitionName, levels = c( "control.control", "control.case" ) ) ) |>
            tidyr::drop_na( TransitionName )
        } )
      }

      if( !is.null( markers ) ) {
        m2td <- m2td[ intersect( markers, names( m2td ) ) ]
      }

      m2td
    },

    # ----- plots -----
    eventsPlot = function(
      showVaxination = TRUE, showCollection = TRUE, showTransitions = TRUE,
      showCovidVariant = TRUE, showSource = TRUE, arrowDaysOffs = 7, arrowCurvature = -0.05,
      subjects = NULL
    ) {
      mm <- self$getMeasures() 
      vaxData <- self$getVaccinationsData()
      mesData <- mm$measurementData |>
        dplyr::distinct( Date, SubjectId, EventCode, isVaxBT, isReinf, VaxReinfCode, Source )
      transData <- mm$transitionData |>
        dplyr::distinct( prevDate, Date, SubjectId, TransitionName )
      transInfo <- mm$transitionInfo

      if( !is.null( subjects ) ) {
        vaxData <- vaxData |> dplyr::mutate( SubjectId = factor( SubjectId, levels = subjects ) ) |> tidyr::drop_na( SubjectId )
        mesData <- mesData |> dplyr::mutate( SubjectId = factor( SubjectId, levels = subjects ) ) |> tidyr::drop_na( SubjectId )
        transData <- transData |> dplyr::mutate( SubjectId = factor( SubjectId, levels = subjects ) ) |> tidyr::drop_na( SubjectId )
      }

      d <- dplyr::bind_rows(
        mesData |> dplyr::select( Date, SubjectId ),
        vaxData |> dplyr::select( Date, SubjectId )
      )

      p <- ggplot2::ggplot( d ) +
        ggplot2::aes( x = Date, y = SubjectId )
      if( showCovidVariant ) {
        p <- p +
          ggplot2::geom_vline( xintercept = as.Date( "2021-06-30" ), linetype = 2, color = "gray" ) +
          ggplot2::geom_vline( xintercept = as.Date( "2022-01-01" ), linetype = 2, color = "gray" )
      }
      if( showVaxination ) {
        p <- p +
          ggplot2::geom_point( size = 1.5, shape = 8, alpha = 0.25, color = "blue", data = vaxData )
      }
      if( showCollection ) {
        p <- p +
          ggplot2::geom_point( 
            data = mesData |> dplyr::filter( isVaxBT | isReinf ) |> dplyr::distinct( Date, SubjectId, VaxReinfCode ), 
            ggplot2::aes( fill = VaxReinfCode ), size = 5, shape = 21, color = "white", alpha = 0.5
          ) +
          ggplot2::scale_fill_manual( 
            values = mm$vaxReinfInfo |> dplyr::pull( VaxReinfColor, VaxReinfCode ) 
          )
      }
      if( showSource ) {
        p <- p +
          ggplot2::geom_point( 
            size = 5, shape = 0, alpha = 1, color = "black", 
            data = mesData |> dplyr::filter( Source == "imm" )
          )
      }
      if( showTransitions ) {
        p <- p +
          ggplot2::geom_curve(
            ggplot2::aes( 
              x = prevDate + arrowDaysOffs, y = SubjectId, 
              xend = Date - arrowDaysOffs, yend = SubjectId, 
              color = TransitionName
            ),
            arrow = ggplot2::arrow( length = ggplot2::unit( 0.125, "cm" ) ),
            data = transData, alpha = 0.5, linewidth = .3,
            curvature = arrowCurvature
          ) +
          ggplot2::scale_color_manual( 
            values = transInfo |> 
              dplyr::distinct( TransitionColor, TransitionName ) |> 
              dplyr::pull( TransitionColor, TransitionName ) 
          )
      }
      if( showCollection ) {
        p <- p +
          ggplot2::geom_text( size = 3, alpha = 1, color = "black", data = mesData, ggplot2::aes( label = EventCode ) )
      }
      p <- p +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid.major=ggplot2::element_line(colour="gray97"),
          panel.grid.minor=ggplot2::element_line(colour="gray99")
        ) +
        ggplot2::scale_x_date( position = "top" )
      p
    },
    markersTimePlots = function( mode, showTransitions = TRUE, showMeasurements = TRUE, transitionAlpha = 0.25, markers = NULL ) {
      mm <- self$getMeasures()
      td <- mm$transitionData
      md <- mm$measurementData

      if( !is.null( markers ) ) {
        td <- td |> dplyr::mutate( Marker = factor( Marker, levels = markers ) ) |> tidyr::drop_na( Marker )
        md <- md |> dplyr::mutate( Marker = factor( Marker, levels = markers ) ) |> tidyr::drop_na( Marker )
      }
      p <- ggplot2::ggplot( td ) +
        ggplot2::aes( x = Date, y = Value )
      if( showTransitions ) {
        p <- p +
          ggplot2::geom_segment(
            ggplot2::aes( x = prevDate, y = prevValue, xend = Date, yend = Value, color = TransitionName ),
            alpha = transitionAlpha, linewidth = .5
          ) +
          ggplot2::scale_color_manual( 
            values = mm$transitionInfo |> dplyr::pull( TransitionColor, TransitionName ) 
          )
      }
      if( is.null( showMeasurements ) ) {
        p <- p +
          ggplot2::geom_point( 
            data = md,
            size = .1, color = "gray30", alpha = .3
          )
      } else if( showMeasurements == "vaxReInf" ) {
        p <- p +
          ggplot2::geom_point( 
            data = md,
            size = .1, color = "gray30", alpha = .3
          ) +
          ggplot2::geom_point( 
            data = md |> dplyr::filter( isVaxBT | isReinf ), 
            ggplot2::aes( fill = VaxReinfCode ), size = 1.5, shape = 21, color = "white", alpha = 0.7
          ) +
          ggplot2::scale_fill_manual( 
            values = mm$vaxReinfInfo |> dplyr::pull( VaxReinfColor, VaxReinfCode ) 
          )
      } else if( showMeasurements == "events" ) {
        p <- p +
          ggplot2::geom_point( 
            data = md, 
            ggplot2::aes( fill = Event ), size = 1.5, shape = 21, color = "white", alpha = 0.7
          ) +
          ggplot2::scale_fill_manual( 
            values = mm$measurementInfo |> dplyr::pull( EventColor, Event ) 
          )
      } else if( showMeasurements == "T1events" ) {
        p <- p +
          ggplot2::geom_point( 
            data = td, 
            ggplot2::aes( fill = prevEvent, y = prevValue, x = prevDate ), size = 1.5, shape = 21, color = "white", alpha = 0.7
          ) +
          ggplot2::scale_fill_manual( 
            values = mm$measurementInfo |> dplyr::pull( EventColor, Event ) 
          )
      } else {
        stop( "Unknown mode (allowed: vaxReInf, events): ", mode)
      }
      p <- p +
        ggplot2::facet_wrap( ~Marker, scales = "free_y" ) +
        ggplot2::theme_bw()
      p
    },
    markersT1T2ScatterPlots = function( markers = NULL ) {
      mm <- self$getMeasures() 
      transData <- mm$transitionData
      transInfo <- mm$transitionInfo

      if( !is.null( markers ) ) {
        transData <- transData |> dplyr::mutate( Marker = factor( Marker, levels = markers ) ) |> tidyr::drop_na( Marker )
      }

      rangeD <- transData |>
        dplyr::group_by( Marker ) |>
        dplyr::summarize( l = min( prevValue, Value ), h = max( prevValue, Value ), .groups = "drop" )
      rangeD <- dplyr::bind_rows(
        rangeD |> dplyr::select( Marker, prevValue = l, Value = l ),
        rangeD |> dplyr::select( Marker, prevValue = h, Value = h )
      )
      p <- ggplot2::ggplot( rangeD ) +
        ggplot2::aes( x = prevValue, y = Value )
      p <- p +
        ggplot2::geom_point( alpha = 0, size = 0.001, color = "blue" ) +
        ggplot2::geom_abline( slope = 1, intercept = 0, linetype = "dashed", color = "gray40" ) +
        ggplot2::geom_point( 
          alpha = 0.75, size = .5, ggplot2::aes( color = TransitionName ), 
          data = transData |> dplyr::sample_n( dplyr::n() )
        ) +
        ggplot2::facet_wrap( ~Marker, scales = "free" ) +
        ggplot2::theme_bw() +
        ggplot2::scale_color_manual( 
          values = transInfo |> 
            dplyr::distinct( TransitionColor, TransitionName ) |> 
            dplyr::pull( TransitionColor, TransitionName ) 
        )
      p
    },
    markersT1T2Stats = function() {
      mm <- self$getMeasures() 
      transData <- mm$transitionData
      transInfo <- mm$transitionInfo

      d <- transData |>
        dplyr::select( Marker, SampleId, Value, prevValue, TransitionName ) |>
        dplyr::left_join(
          mm$measurementData |>
            dplyr::select( Marker, SampleId, isVaxBT ),
          by = c( "Marker", "SampleId" )
        ) |>
        dplyr::mutate( TransitionName = factor( TransitionName, levels = c( "control.case", "control.control" ) ) ) |>
        tidyr::drop_na( TransitionName )

      l <- split( d, d$Marker )
      lapply( l, function( ll ) {
        fit <- glm(TransitionName ~ prevValue, data = ll, family = binomial(link="logit"))
        summary( fit )
      } )
    },
    markersT1forTransitionsBoxPlots = function( markers = NULL ) {
      mm <- self$getMeasures()
      td <- mm$transitionData
      transInfo <- mm$transitionInfo
      md <- mm$measurementData

      if( !is.null( markers ) ) {
        td <- td |> dplyr::mutate( Marker = factor( Marker, levels = markers ) ) |> tidyr::drop_na( Marker )
        md <- md |> dplyr::mutate( Marker = factor( Marker, levels = markers ) ) |> tidyr::drop_na( Marker )
      }
      d <- td |>
        dplyr::select( Marker, SampleId, Value, prevValue, TransitionName ) |>
        tidyr::pivot_longer( c( "Value", "prevValue" ), names_to = "TimePoint", values_to = "Value" ) |>
        dplyr::left_join(
          md |>
            dplyr::select( Marker, SampleId, VaxReinfCode ),
          by = c( "Marker", "SampleId" )
        ) |>
        dplyr::filter( TimePoint == "prevValue" )
      p <- ggplot2::ggplot( d ) +
        ggplot2::aes( x = TransitionName, y = Value, fill = TransitionName )
      p <- p +
        ggplot2::geom_boxplot( outlier.colour = NA, linewidth = 0.5, ggplot2::aes( fill = TransitionName ), color = "black", alpha = 0.1 ) +
        ggplot2::geom_jitter( 
          alpha = 0.25, size = 0.5, color = "black",
          position = ggplot2::position_jitterdodge( 
            jitter.width = 0.3, jitter.height = 0, seed = NA 
          )
        ) +
        ggplot2::facet_wrap( ~ Marker, scales = "free_y" ) +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual( 
          values = transInfo |> 
            dplyr::distinct( TransitionColor, TransitionName ) |> 
            dplyr::pull( TransitionColor, TransitionName ) 
        )
      p
    },
    markersT1T2BoxPlots = function() {
      mm <- self$getMeasures() 
      transData <- mm$transitionData
      transInfo <- mm$transitionInfo

      d <- transData |>
        dplyr::select( Marker, SampleId, Value, prevValue, TransitionName ) |>
        tidyr::pivot_longer( c( "Value", "prevValue" ), names_to = "TimePoint", values_to = "Value" ) |>
        dplyr::left_join(
          mm$measurementData |>
            dplyr::select( Marker, SampleId, VaxReinfCode ),
          by = c( "Marker", "SampleId" )
        ) |>
        dplyr::mutate( VaxReinfCode = dplyr::if_else( TimePoint == "prevValue", ".", VaxReinfCode ) )
      p <- ggplot2::ggplot( d ) +
        ggplot2::aes( x = TimePoint, y = Value, fill = TransitionName, color = VaxReinfCode )
      p <- p +
        ggplot2::geom_boxplot( outlier.colour = NA, linewidth = 0.5, ggplot2::aes( fill = TransitionName ), color = "black", alpha = 0.1 ) +
        ggplot2::geom_jitter( 
          alpha = 0.75, size = 0.5,
          position = ggplot2::position_jitterdodge( 
            jitter.width = 0.1, jitter.height = 0, seed = NA 
          )
        ) +
        ggplot2::facet_wrap( ~ Marker, scales = "free_y" ) +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual( 
          values = transInfo |> 
            dplyr::distinct( TransitionColor, TransitionName ) |> 
            dplyr::pull( TransitionColor, TransitionName ) 
        ) +
        ggplot2::scale_color_manual( 
          values = mm$vaxReinfInfo |> dplyr::pull( VaxReinfColor, VaxReinfCode ) 
        )
      p
    },
    markersBoxPlots = function( markers = NULL ) {
      mm <- self$getMeasures() 
      md <- mm$measurementData

      if( !is.null( markers ) ) {
        md <- md |> dplyr::mutate( Marker = factor( Marker, levels = markers ) ) |> tidyr::drop_na( Marker )
      }

      p <- ggplot2::ggplot( md ) +
        ggplot2::aes( x = Event, y = Value, fill = Event, color = VaxReinfCode ) +
        ggplot2::geom_boxplot( outlier.colour = NA, varwidth = FALSE, color = "black", linewidth = .3, alpha = 0.2 ) +
        ggplot2::geom_jitter( width = 0.2, alpha = 0.75, size = .75 ) +
        ggplot2::facet_wrap( ~Marker, scales = "free_y" ) +
        ggplot2::theme_bw() +
        ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 30, hjust = 1, vjust = 1 ) ) +
        ggplot2::scale_fill_manual( 
          values = mm$measurementInfo |> dplyr::pull( EventColor, Event ) 
        ) +
        ggplot2::scale_color_manual( 
          values = mm$vaxReinfInfo |> dplyr::pull( VaxReinfColor, VaxReinfCode ) 
        )
      p
    },
    vaxReinfT1BoxPlots = function( mode = "vax", splitSymptAsympt = FALSE ) {
      mm <- self$getMeasures() 
      md <- mm$measurementData
      td <- mm$transitionData

      if( mode == "vax" ) {
        d <- md |> dplyr::distinct( SampleId, isVaxBT )
        d <- td |> dplyr::left_join( d, by = "SampleId" ) |> tidyr::drop_na( isVaxBT )
        p <- ggplot2::ggplot( d ) + ggplot2::aes( x = isVaxBT )
      } else if( mode == "reinf" ) {
        d <- md |> dplyr::distinct( SampleId, isReinf )
        d <- td |> dplyr::left_join( d, by = "SampleId" ) |> tidyr::drop_na( isReinf )
        p <- ggplot2::ggplot( d ) + ggplot2::aes( x = isReinf )
      } else {
        stop( "Unknown mode (allowed: vax, reinf): ", mode)
      }
      if( nrow( d ) > 0 ) {
        p <- p +
          ggplot2::aes( y = prevValue ) +
          ggplot2::geom_boxplot( outlier.colour = NA, varwidth = FALSE ) +
          ggplot2::geom_jitter( 
            alpha = 0.75,
            position = ggplot2::position_jitterdodge( jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75, seed = NA ),
            ggplot2::aes( color = Event )
          ) +
          ggplot2::facet_wrap( ~Marker, scales = "free_y" ) +
          ggplot2::theme_bw() +
          ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 30, hjust = 1, vjust = 1 ) ) +
          ggplot2::scale_color_manual( 
            values = mm$measurementInfo |> dplyr::pull( EventColor, Event ) 
          ) +
          ggplot2::scale_fill_manual( 
            values = mm$measurementInfo |> dplyr::pull( EventColor, Event ) 
          )
      }
      p
    },
    ageT1DistPlot = function() {
      mm <- self$getMeasures() 
      d <- mm$transitionData |> 
        dplyr::select( SampleId = prevSampleId, SubjectId ) |> 
        dplyr::left_join( 
          mm$measurementData |> 
          dplyr::distinct( SampleId, AgeYears, Gender ),
          by = "SampleId"
        ) |>
        dplyr::distinct( SubjectId, AgeYears, Gender )
      dd <- d |> 
        dplyr::group_by( Gender ) |>
        dplyr::summarize( n = dplyr::n(), meanAgeYears = mean( AgeYears ), .groups = "drop" )
      p <- ggplot2::ggplot( d ) +
        ggplot2::aes( x = AgeYears ) +
        ggplot2::geom_vline( ggplot2::aes( xintercept = meanAgeYears ), linetype = 2, color = "red", data = dd ) +
        ggplot2::geom_histogram( alpha = 0.5, binwidth = 10, boundary = 50, color = "black" ) +
        ggplot2::geom_rug( alpha = 0.5 ) +
        ggplot2::facet_grid( Gender ~ . ) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle( paste( dd$Gender, "=", dd$n, sep = "", collapse = ", " ) )
      p
    },

    # ----- stats -----
    calcTransitionDataKruskalTest = function( markers = NULL, maxAdjPVal = 0.05 ) {
      m2td <- self$getMarker2transitionData( markers = markers )
      l <- list()
      l$kruskalTests <- lapply( m2td, function( td ) {
        td <- td |> 
          dplyr::mutate( EventCode = factor( EventCode, levels = c( "B", "C", "A", "S" ) ) ) |>
          tidyr::drop_na( EventCode )
        ll <- list()
        ll$td <- td
        tryCatch(
          {
            ll$kTest <- kruskal.test( formula = prevValue ~ EventCode, data = td )
            ll$kTestResult <- do.call( tibble::tibble, ll$kTest )
            ll$pwTest <- pairwise.wilcox.test( td$Value, td$EventCode, p.adjust.method = "BH", exact = FALSE )
            ll$pwTestResult <- ll$pwTest$p.value |> 
              dplyr::as_tibble( rownames = "Term1" ) |> tidyr::gather( "Term2", "adjPVal", -Term1 ) |> tidyr::drop_na()
          },
          error = function( e ) NULL
        )
        ll
      } )

      l$kruskalTestsD <- lapply( l$kruskalTests, function( ll ) ll$kTestResult ) |>
        dplyr::bind_rows( .id = "Marker" ) |>
        dplyr::mutate( adjPVal = p.adjust( p.value, method = "BH" ) ) |>
        dplyr::arrange( adjPVal, p.value ) |>
        dplyr::mutate( sgn = dplyr::if_else( adjPVal <= maxAdjPVal, "*", "" ) ) |>
        dplyr::select( Marker, adjPVal, pVal = p.value, sgn )

      sgnMarkers <- l$kruskalTestsD |> 
        dplyr::filter( adjPVal <= maxAdjPVal ) |> 
        dplyr::pull( Marker )
      
      l$wilcoxTestsD <- lapply( l$kruskalTests[ sgnMarkers ], function( ll ) {
          ll$pwTestResult |>
            dplyr::mutate( sgn = dplyr::if_else( adjPVal <= maxAdjPVal, "*", "" ) )
        } ) |>
        dplyr::bind_rows( .id = "Marker" )
        
      l
    },
    calcTransitionDataWilcoxTest = function( formula, markers = NULL, maxAdjPVal = 0.05 ) {
      m2td <- self$getMarker2transitionData( markers = markers )
      d <- lapply( m2td, function( td ) {
        tryCatch( 
          {
            hTest <- wilcox.test( formula = formula, data = td, exact = FALSE )
            do.call( dplyr::tibble, hTest ) |>
              dplyr::bind_cols( td |> dplyr::count( TransitionName ) |> tidyr::spread( TransitionName, n ) )
          },
          error = function( e ) tibble::tibble()
        )
      } ) |> 
        dplyr::bind_rows( .id = "Marker" ) |>
        dplyr::bind_rows( tibble::tibble( Marker = character(), p.value = numeric(), control.control = integer(), control.case = integer() ) )
      d <- d |>
        dplyr::mutate( adjPVal = p.adjust( p.value, method = "BH" ) ) |>
        dplyr::arrange( adjPVal, p.value ) |>
        dplyr::mutate( sgn = dplyr::if_else( adjPVal <= maxAdjPVal, "*", "" ) ) |>
        dplyr::select( Marker, adjPVal, sgn, pVal = p.value, control.control, control.case )
      d
    },
    calcMeasurementsDataWilcoxTest = function( formula, markers = NULL, maxAdjPVal = 0.05 ) {
      m2md <- self$getMarker2measurementData( markers = markers )
      d <- lapply( m2md, function( md ) {
        md <- md |>
          dplyr::mutate( State = dplyr::if_else( EventCode %in% c( "B", "C" ), "control", "case" ) )
        hTest <- wilcox.test( formula = formula, data = md )
        do.call( dplyr::tibble, hTest ) |>
          dplyr::bind_cols( md |> dplyr::count( State ) |> tidyr::spread( State, n ) )
      } ) |> 
        dplyr::bind_rows( .id = "Marker" ) |>
        dplyr::bind_rows( tibble::tibble( Marker = character(), p.value = numeric(), control = integer(), case = integer() ) )
      d <- d |>
        dplyr::mutate( adjPVal = p.adjust( p.value, method = "BH" ) ) |>
        dplyr::arrange( adjPVal, p.value ) |>
        dplyr::mutate( sgn = dplyr::if_else( adjPVal <= maxAdjPVal, "*", "" ) ) |>
        dplyr::select( Marker, adjPVal, sgn, pVal = p.value, control, case )
      d
    }
  )
)

if( FALSE ) {
  covData <- CovidData$new( 
    "Dataset/BREAK COVID overview v3_combined_data__20240102_LUMC_Mo_SzMK.xlsx", 
    vaxValidDaysRange = c( 14, 270 ), reinfDaysRange = c( 14, 270 ) 
  )
  immData <- ImmunData$new( "Dataset/IMMUN_COVID20240130_V2.xlsx" )
  dm <- DataMerger$new( immData, covData, transDaysRange = c( 14, 14+3*31 ) )
  mm <- dm$getMeasures(); mm
  dp <- DataProcessor$new( dm )
  dp$eventsPlot()
  dp$markersTimePlots( "vaxReInf" )
  dp$markersTimePlots( "events" )
  dp$markersBoxPlots()
  dp$vaxReinfBoxPlots( "vax" )
  dp$vaxReinfBoxPlots( "reinf" )
  dp$markersT1T2BoxPlots()
  dp$markersT1T2ScatterPlots()
}

if( FALSE ) {
  ggplot2::ggplot( mm$measurementData ) +
    ggplot2::aes( x = Date, y = Value, color = Event ) +
    ggplot2::geom_jitter( alpha = 0.9, size = .4, width = 10 ) +
    ggplot2::facet_wrap( ~Marker, scales = "free_y" ) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual( 
      values = mm$measurementInfo |> dplyr::pull( EventColor, Event ) 
    )
}

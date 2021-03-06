% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oasisad_model.R
\name{oasisad_model}
\alias{oasisad_model}
\title{OASISAD model}
\usage{
oasisad_model(
  train_df,
  test_df,
  valid_df = NULL,
  M1 = FALSE,
  refine = FALSE,
  neighbor = FALSE,
  wm_label = NULL,
  re_value = NULL
)
}
\arguments{
\item{train_df}{A data list from oasisad_df function which inlcudes training samples informatin.
If neighbor refinement function will be used, the list should include segmentation and
white matter probability map for each training subject}

\item{test_df}{A data list from oasisad_df function which inlcudes testing samples informatin.
If neighbor refinement function will be used, the list should include segmentation and
white matter probability map for each training subject}

\item{valid_df}{A data list from oasisad_df function which inlcudes validation samples informatin.
If neighbor refinement function will be used, the list should include segmentation and
white matter probability map for each training subject. If it is NULL, optimal threshold
algorithm will be used to calculate threshold}

\item{M1}{A boolean indicates using full model 'M1' or reduced model 'M2',
default is reduced model}

\item{refine}{A boolean incicates whether use OASISAD refinement function,
to refine probability map from logistic regression model}

\item{neighbor}{A boolean incicates whether use neighbor refinement function,
to refine probability map from logistic regression model. If true, segmentation information and
white matter probability of brain is needed}

\item{wm_label}{White matter label in segmentation input}

\item{re_value}{A numeric value will be used in neighor refinement functoin to
refine a voxel's probability of being White matter hyperintensity}
}
\value{
OASISAD model results
}
\description{
The input should be OASISAD data list,
the function will train the model with training and vailidaion data,
then use the testing data to evaluatoin performance
}

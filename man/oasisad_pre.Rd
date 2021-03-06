% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oasisad_pre.R
\name{oasisad_pre}
\alias{oasisad_pre}
\title{OASISAD image preprocessing function}
\usage{
oasisad_pre(
  flair,
  t1,
  t2 = NULL,
  pd = NULL,
  img_space = NULL,
  brain_mask = NULL,
  segmentation = FALSE,
  dir = NULL,
  cores = 1,
  verbose = TRUE
)
}
\arguments{
\item{flair}{Input FLAIR image}

\item{t1}{Input T1 image}

\item{t2}{Input T2 image}

\item{pd}{Input PD image}

\item{img_space}{An image to register and use for creating brain mask if needed.
If NULL, 'flair' image will be used in registration.}

\item{brain_mask}{Input brain_mask, if null, a mask will be obtained by FSL}

\item{segmentation}{A boolean indicates whether using bias_correct for segmentation}

\item{dir}{A user defined output path for fslr segmentation}

\item{cores}{A number indicates how many cores used mclapply}

\item{verbose}{A boolean indicated whether output messages}
}
\value{
preprocessed image data
}
\description{
MRI image data preprocessing with multiple inputs
}
\examples{
library(neurobase)
dl_file = function(url) {
   tfile = file.path(tempdir(), basename(url))
   if (!file.exists(tfile)) {
      if (requireNamespace("httr", quietly = TRUE)) {
          req <- httr::GET(url,
          httr::write_disk(path = tfile))
          httr::stop_for_status(req)
      } else {
         download.file(url, destfile = tfile)
      }
   }
   tfile
}
in_ci <- function() {
 nzchar(Sys.getenv("CI"))
}
on_cran = function() {
 identical(Sys.getenv("NOT_CRAN"), "false")
}
if (in_ci() || on_cran()) {
  if (fslr::have.fsl()) {
    mods = c("FLAIR", "T1W", "T2W", "consensus_gt", "brainmask")
    base_url = file.path(
      "https://raw.githubusercontent.com/muschellij2/open_ms_data",
      "master/cross_sectional/coregistered/patient01/")
    files = paste0(base_url, mods, ".nii.gz")
    files = sapply(files, dl_file)
    names(files) = mods

    flair <- readnii(files["FLAIR"])
    t1 <- readnii(files["T1W"])
    t2 <- readnii(files["T2W"])
    brain_mask <- readnii(files["brainmask"])
    gold_standard = readnii(files["consensus_gt"])
    oasis_preprocessed_data <- oasisad_pre(flair, t1, t2,
      brain_mask = brain_mask, verbose = 2)
  }
}
}

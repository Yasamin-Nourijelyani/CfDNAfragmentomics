library(CfDNAfragmentomics)


test_that("t_test correct false output", {
  controls_bed_file <- system.file("extdata", "d1.bed", package = "CfDNAfragmentomics")
  controls_bed <- read.table(controls_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  sample_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
  sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed)
  ratio

  expect_equal(ratio$cancerous, FALSE)
})


test_that("t_test correct true output", {
  controls_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
  controls_bed <- read.table(controls_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  sample_bed_file <- system.file("extdata", "d1.bed", package = "CfDNAfragmentomics")
  sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed, 0.3, 0.3)
  ratio

  expect_equal(ratio$cancerous, FALSE)
})

test_that("t_test correct true output", {
  controls_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
  controls_bed <- read.table(controls_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  sample_bed_file <- system.file("extdata", "d1.bed", package = "CfDNAfragmentomics")
  sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed)
  ratio

  expect_equal(ratio$cancerous, FALSE)
})



test_that("t-test correct outputs", {

  controls_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
  controls_bed <- read.table(controls_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  sample_bed_file <- system.file("extdata", "d1.bed", package = "CfDNAfragmentomics")
  sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed)
  ratio


  expect_type(ratio, "list")
  expect_length(ratio, 3)
})



#-----------------correct error throwing---------------
test_that("nucleosomeRatio error upon invalid user input", {

  #valid data frames
  controls_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
  controls_bed <- read.table(controls_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  sample_bed_file <- system.file("extdata", "d1.bed", package = "CfDNAfragmentomics")
  sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

  #character second column invalid data frame
  first_column <- c("value_1", "value_2", "value3")
  second_column <- c("value_1", "value_2", "value3")
  third_column <- c(2, 3, 4)
  second_char_df <- data.frame(first_column, second_column, third_column)

  #character third column invalid data frame
  first_column <- c("value_1", "value_2", "value3")
  third_column <- c("value_1", "value_2", "value3")
  second_column <- c(2, 3, 4)
  third_char_df <- data.frame(first_column, second_column, third_column)

  #not enough columns invalid data frame
  first_column <- c("value_1", "value_2", "value3")
  second_column <- c(2, 3, 4)
  not_enough_df <- data.frame(first_column, second_column)



  # large p-value for mono-nucleosome
  expect_error(

    ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed, 5, 0.3)
    )

  # p-value provided as character
  expect_error(
    ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed, 0.5, "0.3")

    )

  # dataframe provided invalid
  expect_error(
    ratio <- nucleosomeRatio(controls_bed = second_char_df, sample_bed = sample_bed, 5, 0.3)


    )

  # dinucleosome provided as character
  expect_error(
    ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed, "0.2", 0.3)
)

  # sample bed provided as character
  expect_error(
    ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = "sample_bed", "0.2", 0.3)

    )

  # sample bed provided invalid df
  expect_error(
    ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed =third_char_df, 0.2, 0.3)

  )

  # sample bed provided invalid df
  expect_error(
    ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed =not_enough_df, 0.2, 0.3)

  )


})

# [END]

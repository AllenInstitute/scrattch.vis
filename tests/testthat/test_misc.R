library(scrattch.vis)

context("misc.R")

test_that("check_genes() Checks input genes against a vector of gene names", {
  library(tasic2016data)
  genes <-
    as.vector(c(
      rownames(tasic2016data::tasic_2016_counts)[1:10],
      "anotherGene",
      "AbADGene-Rik",
      "ViralGene1-20"
    ))
  reference <-
    as.vector(sample(rownames(tasic2016data::tasic_2016_rpkm)))
  output <- check_genes(genes, reference, result = "both")
  expect_true(length(reference) > length(genes))
  expect_true(is.list(output))
  expect_true(length(output) != 0)
  output <- check_genes(genes, reference, result = "matched")
  expect_true(is.character(output))
  expect_false(class(genes) != "character")
  expect_false(class(reference) != "character")
  expect_error(expect_output(check_genes(), NULL))
})



test_that(
  "chr_to_num() Evaluates a character string specifying integer values to a numeric vector",
  {
    test <- "1:5"
    output <- chr_to_num(test)
    expect_true(length(test) == 1)
    expect_true(is.numeric(output))
    expect_error(expect_output(chr_to_num(), NULL))
  }
)


test_that("color_sum() Mixes two colors additively in RGB space",
          {
            expect_true(is.character(color_sum("#1B9E77", "#D95F02")))
            expect_error(color_sum("blue", "white"))
            expect_error(expect_output(color_sum(), NULL))
          })



test_that("pt2mm() Converts font sizes in pt to mm",
          {
            expect_that(pt2mm(12), equals(12 / 2.834645669))
            expect_error(pt2mm(NA))
            expect_error(expect_output(pt2mm(), NULL))
          })



test_that("riken_case() Converts the case of Riken genes no matter input case",
          {
            expect_that(riken_case(c(
              "6330527o06RiK", "A930038C07RIK", "a330070k13rik"
            )), equals(c(
              "6330527O06Rik", "A930038C07Rik", "A330070K13Rik"
            )))
            expect_warning(riken_case("geneA"))
            expect_error(expect_output(riken_case(), NULL))
          })


test_that(
  "split_text() Splits a character string by commas, spaces, tabs, and line breaks into a character vector",
  {
    test <- "Hspa8, Scnn1a,Rbp4    Ptgs2; GeneA:GeneB"
    out <- split_text(test)
    expect_type(test, "character")
    expect_length(test, 1)
    expect_length(out, 6)
    expect_warning(split_text(c("GeneA", "GeneBb2a", "GeneC:2b")))
    expect_error(expect_output(split_text(), NULL))
    
  }
)


test_that("title_case() Converts the case of objects in a character vector to Title Case",
          {
            expect_length(title_case(c("hspa8", "scnn1a", "fhqwghads")), 3)
            expect_type(title_case(c("hspa8", "scnn1a", "fhqwghads")), "character")
            tmp <- tempfile()
            expect_known_output(title_case(c("Scnn1a", "Hspa8")), tmp, print = TRUE)
            expect_warning(title_case(c("hspa8, scnn1a, myrandomgene")))
            expect_error(expect_output(title_case(), NULL))
          })


test_that("values_to_colors() Convert values to colors along a color ramp",
          {
            skip("Pending tests.")
          })
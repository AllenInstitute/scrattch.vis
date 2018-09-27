library(scrattch.vis)

context("Filtering")

test_that("Filter a data.frame with parameters provided in a list", {
  library(tasic2016data)
  anno <- tasic_2016_anno
  filter_list <- list(pass_qc_checks = list(values = "Y",
                                             match_type = "exact"),
                       primary_type = list(values = c("Pvalb","Vip"),
                                           match_type = "grep"))
   
  filtered_anno <- anno %>%
     filter_using_list(filter_list)
  
  expect_true(nrow(anno) == 1809)
  expect_true(nrow(filtered_anno) < nrow(anno))
  expect_true(nrow(filtered_anno) != 0)
  expect_true(class(filter_list) == "list")
  expect_output(str(filter_list) == "match_type")
  expect_output(str(filtered_anno), "pass_qc_checks")
  expect_output(str(filtered_anno), "primary_type")
  expect_that(length(filter_names) == 0, throws_error())
  expect_that(length(missing_cols) > 0, throws_error())
  
})



# Tidyeval notes
library(dplyr)
library(rlang)
library(tasic2016data)

# one character element:
test <- "primary_type_id"

quo_test <- parse_expr(test)

tasic_2016_anno %>%
  group_by(!!quo_test) %>%
  summarise(x = n())
# or
tasic_2016_anno %>%
  group_by(!!parse_expr(test)) %>%
  summarise(x = n())

# one or more character elements
test <- c("primary_type_id","secondary_type_id")

quo_test <- lapply(test, parse_expr)

tasic_2016_anno %>%
  group_by(!!!quo_test) %>%
  summarise(x = n())

# a single entry from multiple elements
test <- c("primary_type_id","secondary_type_id")

quo_test <- lapply(test, parse_expr)

tasic_2016_anno %>%
  group_by(!!quo_test[[1]]) %>%
  summarise(x = n())

# filtering
test <- "primary_type_id"
test_targets <- c(1:5)
tasic_2016_anno %>%
  filter(!!parse_expr(test) %in% test_targets) %>%
  group_by(!!parse_expr(test)) %>%
  summarise(x = n())

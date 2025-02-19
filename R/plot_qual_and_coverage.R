library(tidyvers++e)

import_test_3_w_quality  <- calculate_quality_statistics(import_test_3)

ggplot(import_test_3_w_quality, aes(x = position, y = coverage_no_zero_qual)) + geom_area()
ggplot(import_test_3_w_quality, aes(x = position, y = mean_quality)) + geom_area()

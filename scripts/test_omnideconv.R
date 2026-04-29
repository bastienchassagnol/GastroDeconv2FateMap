library("omnideconv")

data("single_cell_data_1")
data("cell_type_annotations_1")
data("batch_ids_1")
data("bulk")

single_cell_data <- single_cell_data_1[1:2000, 1:500]
cell_type_annotations <- cell_type_annotations_1[1:500]
batch_ids <- batch_ids_1[1:500]
bulk <- bulk[1:2000, ]


deconv_bisque <- deconvolute(
  bulk,
  NULL,
  "bisque",
  single_cell_data,
  cell_type_annotations,
  batch_ids
)

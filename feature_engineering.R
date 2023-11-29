# the L1 distances between axonal and dendritic in x/y/z
# mention that i know these distances are imperfect because of the footnote in the project description doc
# say i don't know much about this neuron stuff but i'm assuming abs is the way to go
L1_ad_xyz <- function(df) {
  df |>
    dplyr::mutate(ad_dx = abs(axonal_coor_x - dendritic_coor_x),
                  ad_dy = abs(axonal_coor_y - dendritic_coor_y),
                  ad_dz = abs(axonal_coor_z - dendritic_coor_z))
}

# axonal to pre_nucleus and dendritic to post_nucleus distances
ad_to_nucleus <- function(df) {
  df |>
    dplyr::mutate(apre_dx = abs(axonal_coor_x - pre_nucleus_x),
                  apre_dy = abs(axonal_coor_y - pre_nucleus_y),
                  apre_dz = abs(axonal_coor_z - pre_nucleus_z),
                  apre_d = sqrt(apre_dx^2 + apre_dy^2 + apre_dz^2),
                  dpost_dx = abs(dendritic_coor_x - post_nucleus_x),
                  dpost_dy = abs(dendritic_coor_y - post_nucleus_y),
                  dpost_dz = abs(dendritic_coor_z - post_nucleus_z),
                  dpost_d = sqrt(dpre_dx^2 + dpre_dy^2 + dpre_dz^2))
}

# axonal to post_nucleus and dendritic to pre_nucleus distances (maybe this distance metric is useful?)
ad_to_other_nucleus <- function(df) {
  df |>
    dplyr::mutate(apost_dx = abs(axonal_coor_x - post_nucleus_x),
                  apost_dy = abs(axonal_coor_y - post_nucleus_y),
                  apost_dz = abs(axonal_coor_z - post_nucleus_z),
                  apost_d = sqrt(apre_dx^2 + apre_dy^2 + apre_dz^2),
                  dpre_dx = abs(dendritic_coor_x - pre_nucleus_x),
                  dpre_dy = abs(dendritic_coor_y - pre_nucleus_y),
                  dpre_dz = abs(dendritic_coor_z - pre_nucleus_z),
                  dpre_d = sqrt(dpre_dx^2 + dpre_dy^2 + dpre_dz^2))
}

# pre_nucleus to post_nucleus
pre_to_post_nucleus <- function(df) {
  df |>
    dplyr::mutate(prepost_dx = abs(pre_nucleus_x - post_nucleus_x),
                  prepost_dy = abs(pre_nucleus_y - post_nucleus_y),
                  prepost_dz = abs(pre_nucleus_z - post_nucleus_z),
                  prepost_d = sqrt(prepost_dx^2 + prepost_dy^2 + prepost_dz^2))
}

# cosine similarity for each
morph_feature_sim <- function(df) {
  # gonna need to pivot longer here perhaps
}

# abs(pre_morph_emb_0 - post_morph_emb_0) essentially
morph_diff <- function(df) {
  # gonna need to pivot longer here perhaps
}

# abs(pre_feature_weight_0 - post_feature_weight_0) essentially
feature_diff <- function(df) {
  # gonna need to pivot longer here perhaps
}


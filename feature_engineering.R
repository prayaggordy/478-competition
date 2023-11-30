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
                  dpost_d = sqrt(dpost_dx^2 + dpost_dy^2 + dpost_dz^2))
}

# axonal to post_nucleus and dendritic to pre_nucleus distances (maybe this distance metric is useful?)
ad_to_other_nucleus <- function(df) {
  df |>
    dplyr::mutate(apostnuc_dx = abs(axonal_coor_x - post_nucleus_x),
                  apostnuc_dy = abs(axonal_coor_y - post_nucleus_y),
                  apostnuc_dz = abs(axonal_coor_z - post_nucleus_z),
                  apostnuc_d = sqrt(apostnuc_dx^2 + apostnuc_dy^2 + apostnuc_dz^2),
                  dprenuc_dx = abs(dendritic_coor_x - pre_nucleus_x),
                  dprenuc_dy = abs(dendritic_coor_y - pre_nucleus_y),
                  dprenuc_dz = abs(dendritic_coor_z - pre_nucleus_z),
                  dprenuc_d = sqrt(dprenuc_dx^2 + dprenuc_dy^2 + dprenuc_dz^2))
}

# pre_nucleus to post_nucleus
pre_to_post_nucleus <- function(df) {
  df |>
    dplyr::mutate(prepost_nuc_dx = abs(pre_nucleus_x - post_nucleus_x),
                  prepost_nuc_dy = abs(pre_nucleus_y - post_nucleus_y),
                  prepost_nuc_dz = abs(pre_nucleus_z - post_nucleus_z),
                  prepost_nuc_d = sqrt(prepost_nuc_dx^2 + prepost_nuc_dy^2 + prepost_nuc_dz^2))
}

pre_to_post_rf <- function(df) {
  df |>
    dplyr::mutate(prepost_rf_dx = abs(pre_rf_x - post_rf_x),
                  prepost_rf_dy = abs(pre_rf_y - post_rf_y),
                  prepost_rf_d = sqrt(prepost_rf_dx^2 + prepost_rf_dy^2))
}

projection_region <- function(df) {
  df |>
    dplyr::mutate(brain_area = paste(pre_brain_area, post_brain_area, sep = "_"))
}

# abs(pre_morph_emb_0 - post_morph_emb_0) essentially
# and the cosine similarity
diff_and_sim <- function(df, prefix) {
  pivoted <- df_joined |>
    dplyr::select(pre_nucleus_id, post_nucleus_id,
                  dplyr::contains(prefix)) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(cols = dplyr::contains(prefix),
                        names_to = c(".value", "num"),
                        names_pattern = "(pre|post)_(.*)")

  sim <- pivoted |>
    dplyr::group_by(pre_nucleus_id, post_nucleus_id) |>
    dplyr::summarize("{prefix}_sim" := lsa::cosine(pre, post) |> as.vector(),
                     "pre_{prefix}_compactness" := norm(pre, type = "2"),
                     "post_{prefix}_compactness" := norm(post, type = "2"))

  pivoted |>
    dplyr::mutate(diff = pre - post) |>  # probably not abs here
    dplyr::select(-c(pre, post)) |>
    tidyr::pivot_wider(names_from = num, names_prefix = "diff_",
                       values_from = diff) |>
    dplyr::left_join(sim, by = c("pre_nucleus_id", "post_nucleus_id"))
}

e_dist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2, na.rm = T))
}

fe_distances <- function(df) {
  df |>
    dplyr::select(ID, dplyr::ends_with(c("_x", "_y", "_z"))) |>
    tidyr::pivot_longer(cols = -ID,
                        names_to = c(".value", "coordinate"),
                        names_pattern = "(.*)_(x|y|z)") |>
    dplyr::group_by(ID) |>
    dplyr::summarize(dist_axonal_postN = e_dist(axonal_coor, post_nucleus),
                     dist_dendritic_preN = e_dist(dendritic_coor, pre_nucleus),
                     dist_preN_postN = e_dist(pre_nucleus, post_nucleus),
                     dist_preRF_postRF = e_dist(pre_rf, post_rf))
}

projection_regions <- function(df) {
  df |>
    dplyr::mutate(projection_region = paste(pre_brain_area, post_brain_area, sep = "_"),
                  projection_region_change = pre_brain_area != post_brain_area) |>
    dplyr::select(ID, projection_region, projection_region_change)
}

count_metrics <- function(df) {
  within_pct <- df |>
    dplyr::count(pre_nucleus_id, post_nucleus_id, name = "within_pct") |>
    dplyr::mutate(within_pct = dplyr::percent_rank(within_pct))

  pre_without_pct <- df |>
    dplyr::group_by(pre_nucleus_id) |>
    dplyr::summarize(pre_within_pct = dplyr::n_distinct(post_nucleus_id)) |>
    dplyr::mutate(pre_within_pct = dplyr::percent_rank(pre_within_pct))

  post_without_pct <- df |>
    dplyr::group_by(post_nucleus_id) |>
    dplyr::summarize(post_without_pct = dplyr::n_distinct(pre_nucleus_id)) |>
    dplyr::mutate(post_without_pct = dplyr::percent_rank(post_without_pct))

  within_pct |>
    dplyr::left_join(pre_without_pct, by = "pre_nucleus_id") |>
    dplyr::left_join(post_without_pct, by = "post_nucleus_id")
}

impute_me <- function(df, fe_type, fit_fn = "imputation/impute_fit.Rds") {
  if (fe_type == "training") {
    big_fit <- lm(formula(paste0("cbind(",
                                 paste("pre_morph_emb", 0:31, sep = "_", collapse = ", "),
                                 ") ~ .")),
                  data = df |>
                    dplyr::select(-dplyr::all_of(c("ID", "pre_nucleus_id", "post_nucleus_id", "connected"))))
    save(big_fit, file = fit_fn)
  } else {
    load(fit_fn)
  }

  preds <- predict(big_fit, df) |>
    tibble::as_tibble() |>
    dplyr::mutate(ID = df$ID) |>
    tidyr::pivot_longer(cols = -ID,
                        values_to = "pred")

  df |>
    dplyr::select(ID, dplyr::starts_with("pre_morph_emb")) |>
    tidyr::pivot_longer(cols = -ID) |>
    dplyr::inner_join(preds, by = c("ID", "name")) |>
    dplyr::mutate(v = dplyr::coalesce(value, pred)) |>
    dplyr::select(ID, name, v) |>
    tidyr::pivot_wider(names_from = name, values_from = v)
}

compress_fw_me <- function(df, prefix) {
  df |>
    dplyr::select(pre_nucleus_id, post_nucleus_id,
                  dplyr::contains(prefix)) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(cols = dplyr::contains(prefix),
                        names_to = c(".value", "num"),
                        names_pattern = "(pre|post)_(.*)") |>
    dplyr::mutate(diff = pre - post) |>
    dplyr::group_by(pre_nucleus_id, post_nucleus_id) |>
    dplyr::summarize("{prefix}_sim" := lsa::cosine(pre, post) |> as.vector(),
                     "pre_{prefix}_norm" := norm(pre, type = "2"),
                     "post_{prefix}_norm" := norm(post, type = "2"),
                     "diff_{prefix}_norm" := norm(diff, type = "2")) |>
    dplyr::ungroup()
}

do_feature_engineering <- function(df_raw, fe_type,
                                   df_fw = feature_weights,
                                   df_me = morph_embeddings) {
  first_features <- df_raw |>
    dplyr::relocate(c(pre_nucleus_id, post_nucleus_id), .after = ID) |>
    dplyr::left_join(fe_distances(df_raw), by = "ID") |>
    dplyr::left_join(projection_regions(df_raw), by = "ID") |>
    dplyr::left_join(count_metrics(df_raw), by = c("pre_nucleus_id", "post_nucleus_id"))  # training/testing split will happen first

  morph_imputed <- first_features |>
    dplyr::left_join(df_me |>
                       dplyr::rename_with(.fn = \(x) paste("pre", x, sep = "_")),
                     by = "pre_nucleus_id") |>
    impute_me(fe_type = fe_type)

  df_big <- first_features |>
    dplyr::left_join(morph_imputed, by = "ID") |>
    dplyr::left_join(df_me |>
                       dplyr::rename_with(.fn = \(x) paste("post", x, sep = "_")),
                     by = "post_nucleus_id") |>
    dplyr::left_join(df_fw |>
                       dplyr::rename_with(.fn = \(x) paste("pre", x, sep = "_")),
                     by = "pre_nucleus_id") |>
    dplyr::left_join(df_fw |>
                       dplyr::rename_with(.fn = \(x) paste("post", x, sep = "_")),
                     by = "post_nucleus_id")

  compressed <- df_big |>
    compress_fw_me("morph_emb") |>
    dplyr::left_join(
      df_big |>
        compress_fw_me("feature_weight"),
      by = c("pre_nucleus_id", "post_nucleus_id")
    )

  first_features |>
    dplyr::left_join(compressed, by = c("pre_nucleus_id", "post_nucleus_id"))
}

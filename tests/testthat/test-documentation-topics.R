find_package_root <- function(){
  current <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if(file.exists(file.path(current, "DESCRIPTION"))){
      return(current)
    }
    parent <- dirname(current)
    if(identical(parent, current)){
      stop("Cannot locate package root.")
    }
    current <- parent
  }
}

read_man_file <- function(topic){
  root <- tryCatch(find_package_root(), error = function(e) NULL)
  if(!is.null(root)){
    path <- file.path(root, "man", paste0(topic, ".Rd"))
    if(file.exists(path)){
      return(readLines(path, warn = FALSE))
    }
  }

  rd <- read_installed_rd(topic)
  deparse(rd)
}

read_installed_rd <- function(topic){
  rd <- tools::Rd_db("BayesGP")[[paste0(topic, ".Rd")]]
  if(is.null(rd)){
    stop("Cannot locate Rd topic: ", topic)
  }
  rd
}

read_man_aliases <- function(topic){
  root <- tryCatch(find_package_root(), error = function(e) NULL)
  if(!is.null(root)){
    path <- file.path(root, "man", paste0(topic, ".Rd"))
    if(file.exists(path)){
      rd <- tools::parse_Rd(path)
    } else {
      rd <- read_installed_rd(topic)
    }
  } else {
    rd <- read_installed_rd(topic)
  }

  alias_nodes <- rd[vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\alias"), logical(1))]
  unlist(lapply(alias_nodes, as.character), use.names = FALSE)
}

man_topic_exists <- function(topic){
  root <- tryCatch(find_package_root(), error = function(e) NULL)
  if(!is.null(root)){
    return(file.exists(file.path(root, "man", paste0(topic, ".Rd"))))
  }
  paste0(topic, ".Rd") %in% names(tools::Rd_db("BayesGP"))
}

test_that("model-specific f documentation aliases are generated", {
  expected_aliases <- list(
    f_iid = c("f_iid", "f.iid"),
    f_iwp = c("f_iwp", "f.iwp"),
    f_sgp = c("f_sgp", "f.sgp"),
    f_mgp = c("f_mgp", "f.mgp"),
    f_tiwp2 = c("f_tiwp2", "f.tiwp2")
  )

  for(topic in names(expected_aliases)){
    aliases <- read_man_aliases(topic)
    for(alias in expected_aliases[[topic]]){
      expect_true(alias %in% aliases, info = alias)
    }
  }
})

test_that("generated f documentation avoids legacy FEM wording", {
  topics <- c("f", "f_iid", "f_iwp", "f_sgp", "f_mgp", "f_tiwp2", "model_fit")
  rd <- unlist(lapply(topics, read_man_file), use.names = FALSE)

  expect_false(any(grepl("legacy FEM", rd, fixed = TRUE)))
  expect_false(any(grepl("legacy O-spline", rd, fixed = TRUE)))
})

test_that("model_arguments is not documented", {
  expect_false(man_topic_exists("model_arguments"))
  expect_false("model_arguments" %in% getNamespaceExports("BayesGP"))
})

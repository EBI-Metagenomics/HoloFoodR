let
  pkgs =
    import
      (fetchTarball "https://github.com/rstats-on-nix/nixpkgs/archive/35097b60acd597bcbb34698b6691dd5be2e39e3c.tar.gz")
      { };

  rpkgs = builtins.attrValues {
    inherit (pkgs.rPackages)
      ALDEx2
      bartMachine
      basilisk
      BiocManager
      BiocVersion
      biomformat
      chromote
      ComplexHeatmap
      devtools
      dplyr
      DT
      FSA
      gbp
      GGally
      ggplot2
      ggpubr
      glmnet
      knitr
      latex2exp
      MGnifyR
      miaViz
      miniUI
      MOFA2
      multiview
      patchwork
      pheatmap
      randomForest
      renv
      reshape2
      reticulate
      rmarkdown
      shadowtext
      styler
      UpSetR
      vegan
      xgboost
      languageserver
      ;
  };

  git_archive_pkgs = [
    (pkgs.rPackages.buildRPackage {
      name = "rix";
      src = pkgs.fetchgit {
        url = "https://github.com/b-rodrigues/rix/";
        rev = "008e3d1ac579cbcd4de3b9a2315635c652b687bd";
        sha256 = "sha256-nrHHiHLYXehly0/qwhGvGmIgJrgyiW4kFIKJ1HvWQfg=";
      };
      propagatedBuildInputs = builtins.attrValues {
        inherit (pkgs.rPackages)
          codetools
          curl
          jsonlite
          sys
          ;
      };
    })

    (pkgs.rPackages.buildRPackage {
      name = "HoloFoodR";
      src = pkgs.fetchgit {
        url = "https://github.com/EBI-Metagenomics/HoloFoodR";
        rev = "93728e327ac50f950a01fdf8af359ccccf70d890";
        sha256 = "sha256-hXWDpICPCZ33jgX3ntAy2wnhyVFptaOGkWHzQaHgQuY=";
      };
      propagatedBuildInputs = builtins.attrValues {
        inherit (pkgs.rPackages)
          TreeSummarizedExperiment
          MultiAssayExperiment
          dplyr
          httr2
          jsonlite
          S4Vectors
          ;
      };
    })

    (pkgs.rPackages.buildRPackage {
      name = "mia";
      src = pkgs.fetchgit {
        url = "https://github.com/microbiome/mia";
        rev = "52d2a236449206f66e8b03ad21e7abe40313c2cf";
        sha256 = "sha256-Y40PB8/rFwite7cvRTbz6qbjj3wEli+M0oFt/qI4Q5E=";
      };
      propagatedBuildInputs = builtins.attrValues {
        inherit (pkgs.rPackages)
          SummarizedExperiment
          SingleCellExperiment
          TreeSummarizedExperiment
          MultiAssayExperiment
          MASS
          ape
          decontam
          vegan
          BiocGenerics
          S4Vectors
          IRanges
          Biostrings
          DECIPHER
          BiocParallel
          DelayedArray
          DelayedMatrixStats
          scuttle
          scater
          DirichletMultinomial
          rlang
          dplyr
          tibble
          tidyr
          bluster
          MatrixGenerics
          mediation
          rbiom
          ;
      };
    })
  ];

  R = with pkgs; [
    (rWrapper.override {
      packages = [
        rpkgs
        git_archive_pkgs
      ];
    })
  ];

  rstudio = with pkgs; [
    (rstudioWrapper.override {
      packages = [
        rpkgs
        git_archive_pkgs
      ];
    })
  ];

  system_packages = builtins.attrValues {
    inherit (pkgs)
      glibcLocales
      nix
      ;
  };

in

pkgs.mkShell {
  LOCALE_ARCHIVE =
    if pkgs.system == "x86_64-linux" then "${pkgs.glibcLocales}/lib/locale/locale-archive" else "";
  LANG = "en_US.UTF-8";
  LC_ALL = "en_US.UTF-8";
  LC_TIME = "en_US.UTF-8";
  LC_MONETARY = "en_US.UTF-8";
  LC_PAPER = "en_US.UTF-8";
  LC_MEASUREMENT = "en_US.UTF-8";

  buildInputs = [
    R
    rstudio
    system_packages
  ];

}

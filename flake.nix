{
  description = "Nix Flake HoloFoodR package";

  inputs = {
    nixpkgs.url = "github:rstats-on-nix/nixpkgs/8d9fcbc47a3574dc897f84fa48fa3cb8f7771761"; # Bioconductor devel
    flake-utils.url = "github:numtide/flake-utils";
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
        HoloFoodR = pkgs.rPackages.buildRPackage {
          name = "HoloFoodR";
          src = self;
          propagatedBuildInputs = builtins.attrValues {
            inherit (pkgs.rPackages)
              dplyr
              httr2
              jsonlite
              S4Vectors
              BiocStyle
              DT
              ggh4x
              ggsignif
              knitr
              MGnifyR
              mia
              miaViz
              MOFA2
              patchwork
              reticulate
              rmarkdown
              scater
              shadowtext
              testthat
              UpSetR
              ;
          };
        };
        R = with pkgs; [
          (rWrapper.override {
            packages = [
              HoloFoodR
            ];
          })
        ];
        system_packages = builtins.attrValues {
          inherit (pkgs)
            quarto
            glibcLocales
            nix
            ;
        };
      in
      {
        devShells.default = pkgs.mkShell {
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
            system_packages
          ];
        };
      }
    );
}

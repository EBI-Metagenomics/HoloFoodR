{
  inputs = {
    nixpkgs.url = "github:rstats-on-nix/nixpkgs/c2a9a8c6b785e80d0188c29f2dc9f1870c95b812";
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
        renv = with pkgs.rPackages; [
          dplyr
          HoloFoodR
        ];
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            R
            renv
          ];
        };
      }
    );
}

{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs";
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

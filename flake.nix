{
  inputs = {
    nixpkgs.url = "github:rstats-on-nix/nixpkgs/c2a9a8c6b785e80d0188c29f2dc9f1870c95b812";

  };
  outputs =
    { self, nixpkgs }:
    let
      supportedSystems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
      forEachSupportedSystem =
        f:
        nixpkgs.lib.genAttrs supportedSystems (
          system:
          f {
            pkgs = import nixpkgs {
              inherit system;
              overlays = [ self.overlays.default ];
            };
          }
        );
    in
    {
      overlays.default = final: prev: rec {
        rEnv = final.rWrapper.override {
          packages = with final.rPackages; [ HoloFoodR ];
        };
      };

      devShells = forEachSupportedSystem (
        { pkgs }:
        {
          default = pkgs.mkShell {
            packages = with pkgs; [
              rEnv
              quarto
            ];
          };
        }
      );
    };
}

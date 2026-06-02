{
  pkgs ?
    import (
      builtins.fetchTarball {
        name = "nixos-gcc13";
        url = "https://github.com/nixos/nixpkgs/archive/e24b4c09e963677b1beea49d411cd315a024ad3a.tar.gz";
        sha256 = "1m383nldy1jvl8n2jw02l26w9iib4b2a9341c1kf74miaahw7qx6";
      }
    )
    {},
}:
pkgs.callPackage ./package.nix {}

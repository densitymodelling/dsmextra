.onAttach <- function(libname, pkgname) {
  packageStartupMessage("-----------------------------------------------\n",
                        "dsmextra: version 1.0.0:\n",
                        "-----------------------------------------------\n",
                        "* Please cite as:\nBouchet PJ, Miller DL, Mannocci L (2019). A Toolkit for Extrapolation Assessments\nin Density Surface Models. R package version 1.0.0.\n",
                        "\n* Quick start guide:\nA vignette is available at:\nhttps://densitymodelling.github.io/model-extrapolation/vignette/Extrapolation-vignette.html\n",
                        "\n* Further information:\ndsmextra is an output of the DenMod project. For more details, please visit:\nhttps://synergy.st-andrews.ac.uk/denmod/")
}

import os

SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")

# ====================================================================
# Include modular rule files
# ====================================================================

include: "bam_refiner.smk"
include: "methylation.smk"
include: "copynumber.smk"
include: "nanomonsv.smk"
include: "clairs.smk"
include: "deepsomatic.smk"

rm containers/pasa/flag_pasa.image
singularity pull flag_pasa.image docker://ghcr.io/formbio/flag_pasa:latest
mv flag_pasa.image containers/pasa/

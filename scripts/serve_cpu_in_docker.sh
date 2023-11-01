if [ -z "${ASKCOS_REGISTRY}" ]; then
  export ASKCOS_REGISTRY=registry.gitlab.com/mlpds_mit/askcosv2
fi

docker run -d --rm \
  --name mcts \
  --network=host \
  -t ${ASKCOS_REGISTRY}/tree_search/mcts:1.0-cpu

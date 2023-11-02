# MCTS (for Tree builder)

---
Serving module for MCTS Tree builder, based on ASKCOSv1 tree builder v2

## Step 1: Environment Setup

First set up the url to the remote registry
```
export ASKCOS_REGISTRY=registry.gitlab.com/mlpds_mit/askcosv2/askcos2_core
```

### Using Docker

- Only option: build from local
```
docker build -f Dockerfile_cpu -t ${ASKCOS_REGISTRY}/tree_search/mcts:1.0-cpu .
```

### Using Singularity

- Only option: build from local
```
singularity build -f mcts_cpu.sif singularity_cpu.def
```

## Step 2: Start the Service

### Using Docker

```
sh scripts/serve_cpu_in_docker.sh
```

### Using Singularity

```
sh scripts/serve_cpu_in_singularity.sh
```

Note that these scripts start the service in the background (i.e., in detached mode). So they would need to be explicitly stopped if no longer in use
```
(Docker)        docker stop mcts
(Singularity)   singularity instance stop mcts
```

## Step 3: Query the Service
Note that the expand_one service depends on the API gateway in askcos2_core (quite a number of services). Therefore, the query would only work if the API gateway (with the backend services) has been started.

Sample Query
```
curl http://0.0.0.0:9301/get_outcomes \
	--header "Content-Type: application/json" \
	--request POST \
	--data '{"original": "CC(=O)c1ccc2c(ccn2C(=O)OC(C)(C)C)c1", "outcomes": ["CC(=O)c1ccc2[nH]ccc2c1.CC(C)(C)OC(=O)OC(=O)OC(C)(C)C", "CC(=O)Cl.CC(C)(C)OC(=O)n1ccc2ccccc21", "CON(C)C(=O)c1ccc2c(ccn2C(=O)OC(C)(C)C)c1.C[Mg+]", "CC(C)(C)OC(=O)n1ccc2cc(Br)ccc21.CON(C)C(C)=O", "CC(O)c1ccc2c(ccn2C(=O)OC(C)(C)C)c1", "CC(O)c1ccc2c(ccn2C(=O)OC(C)(C)C)c1", "CC(=O)OC(C)=O.CC(C)(C)OC(=O)n1ccc2ccccc21", "CC(=O)Cl.CC(C)(C)OC(=O)n1ccc2cc(Br)ccc21"], "cluster_method": "hdbscan"}'

```
Sample response
```
{
    "status":"SUCCESS",
    "error":"",
    "results":[
        [0,1,2,3,4,5,6,7],
        {
            "0":"Reaction Cluster #1",
            "1":"Reaction Cluster #2",
            "2":"Reaction Cluster #3",
            "3":"Reaction Cluster #4",
            "4":"Reaction Cluster #5",
            "5":"Reaction Cluster #6",
            "6":"Reaction Cluster #7",
            "7":"Reaction Cluster #8"
        }
    ]
}
```

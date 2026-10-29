#!/usr/bin/env python3
"""MPI-friendly launcher: assign local ranks round-robin over visible GPUs."""
import argparse
import os

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--gpus-per-node", type=int, required=True)
parser.add_argument("--backend", choices=("auto", "cpu", "cuda"), default="auto")
parser.add_argument("command", nargs=argparse.REMAINDER)
args = parser.parse_args()
if args.gpus_per_node < 1:
    parser.error("--gpus-per-node must be positive")
command = args.command[1:] if args.command[:1] == ["--"] else args.command
if not command:
    parser.error("provide the NWChem executable and input after --")
rank = next((os.environ[key] for key in
             ("OMPI_COMM_WORLD_LOCAL_RANK", "MPI_LOCALRANKID", "SLURM_LOCALID", "PMI_LOCAL_RANK")
             if key in os.environ), "0")
os.environ.setdefault("NWGPU_DEVICE", str(int(rank) % args.gpus_per_node))
os.environ["NWGPU_BACKEND"] = args.backend
os.execvp(command[0], command)

# Testing the Devcontainer Locally

## Prerequisites

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
2. Install [VS Code](https://code.visualstudio.com/)
3. Install the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

## Testing Steps

1. Open this project in VS Code:
   ```bash
   code /path/to/nextflow-curso-herramientas-computacionales-2025
   ```

2. Open the Command Palette (Cmd+Shift+P on Mac, Ctrl+Shift+P on Windows/Linux)

3. Type and select: **"Dev Containers: Rebuild and Reopen in Container"**

4. Wait for the container to build (first time will take several minutes)

5. Once the container is ready, you'll see the setup script output in the terminal

6. Verify the environment:
   ```bash
   nextflow -version
   docker --version
   micromamba --version
   ```

## Troubleshooting

If the build fails:

1. Check Docker is running: `docker ps`
2. View build logs in the VS Code terminal
3. Try: **"Dev Containers: Rebuild Container"** to start fresh

## Testing Changes

After modifying devcontainer files:
- Run **"Dev Containers: Rebuild Container"** to apply changes
- Use **"Dev Containers: Reopen Folder Locally"** to exit the container

## Alternative: Build Without VS Code

You can also build the container manually:

```bash
cd .devcontainer
docker build -t nextflow-curso:test .
docker run -it --rm -v $(pwd)/..:/workspace nextflow-curso:test bash
```
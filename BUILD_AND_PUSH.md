# Building and publishing the Traj2 Docker image

The image is published to GitHub Container Registry (GHCR) at:

```
ghcr.io/ru-aging/traj2:<version>
```

## Automated path (recommended)

The `.github/workflows/docker-publish.yml` workflow builds and pushes
the image automatically when you push a version tag.

```bash
# 1. Make sure Dockerfile, docker/entrypoint.sh, and macros are committed.
git add Dockerfile docker/ docker-compose.yml .github/workflows/docker-publish.yml
git commit -m "Add Docker image build pipeline"
git push origin main

# 2. Tag a release.
git tag v1.0.0
git push origin v1.0.0
```

The workflow takes ~5–10 minutes. Watch progress in the **Actions** tab on
GitHub. When it succeeds, the image appears under
**Packages** on the repo's main page.

## One-time setup

1. **Make the package public** (optional but expected for an open paper):
   - GitHub → repo → Packages → `traj2` → Package settings → Change visibility → Public.

2. **Confirm Actions has package-write permission**:
   - GitHub → repo → Settings → Actions → General → Workflow permissions:
     select "Read and write permissions" and save. (This is required so
     `secrets.GITHUB_TOKEN` can push to GHCR.)

## Manual path (if you need to push from a laptop)

You shouldn't normally need this, but for the record:

```bash
# 1. Authenticate. Use a Personal Access Token with write:packages scope.
echo "$GHCR_PAT" | docker login ghcr.io -u <your-github-username> --password-stdin

# 2. Build, tag, push.
docker build -t ghcr.io/ru-aging/traj2:1.0.0 .
docker push ghcr.io/ru-aging/traj2:1.0.0

# Optional 'latest' tag.
docker tag ghcr.io/ru-aging/traj2:1.0.0 ghcr.io/ru-aging/traj2:latest
docker push ghcr.io/ru-aging/traj2:latest
```

## Verifying the published image

```bash
docker pull ghcr.io/ru-aging/traj2:1.0.0
docker run --rm ghcr.io/ru-aging/traj2:1.0.0 help
```

`validate` and `run` require a SAS installation mounted at `/opt/sas`,
so those commands cannot be smoke-tested without it.

## Troubleshooting

- **403 on push**: package permissions not set. See "One-time setup" #2.
- **Image name rejected by GHCR**: must be lowercase. The workflow handles
  this through `${{ github.repository_owner }}`, which returns the org name
  as-is. If your org has uppercase letters and GHCR rejects them, hardcode
  `IMAGE_NAME: ru-aging/traj2` in `docker-publish.yml`.
- **"unauthorized" pulling the public image**: GHCR sometimes requires
  `docker logout ghcr.io` before pulling a public image without a token.

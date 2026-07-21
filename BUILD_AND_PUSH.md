# Building and publishing the Traj2 Docker image

The image is published to Docker Hub at:

```
chaolab/gbtm-macros:<version>
```

## Prerequisites (one-time)

1. A Docker Hub organization named `chaolab` with a repository
   `gbtm-macros`. Create the repository on Docker Hub and set it to **Public**.
   (Docker Hub does not auto-create org repositories on first push.)
2. Two GitHub Actions secrets on this repo
   (**Settings → Secrets and variables → Actions**):
   - `DOCKERHUB_USERNAME` : a Docker Hub account with push access to the `chaolab` org.
   - `DOCKERHUB_TOKEN` : a Docker Hub **access token** (not the account password).

## Automated path (recommended)

The `.github/workflows/docker-publish.yml` workflow builds and pushes the
image automatically when you push a version tag.

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
GitHub. When it succeeds, the image appears on Docker Hub under
`chaolab/gbtm-macros`.

## Manual path (if you need to push from a laptop)

You shouldn't normally need this, but for the record:

```bash
# 1. Authenticate to Docker Hub. Use a Docker Hub access token.
echo "$DOCKERHUB_TOKEN" | docker login -u <your-dockerhub-username> --password-stdin

# 2. Build, tag, push.
docker build -t chaolab/gbtm-macros:1.0.0 .
docker push chaolab/gbtm-macros:1.0.0

# Optional 'latest' tag.
docker tag chaolab/gbtm-macros:1.0.0 chaolab/gbtm-macros:latest
docker push chaolab/gbtm-macros:latest
```

## Verifying the published image

```bash
docker pull chaolab/gbtm-macros:1.0.0
docker run --rm chaolab/gbtm-macros:1.0.0 help
```

`validate` and `run` require a SAS installation mounted at `/opt/sas`,
so those commands cannot be smoke-tested without it.

## Troubleshooting

- **401/403 on push**: the `DOCKERHUB_USERNAME` / `DOCKERHUB_TOKEN` secrets are
  missing or expired, or the account lacks push access to the `chaolab` org.
- **"repository does not exist"**: create `chaolab/gbtm-macros` on Docker Hub
  first (Prerequisites #1).
- **Image name rejected**: Docker Hub names must be lowercase; `chaolab/gbtm-macros`
  already is.
- **"unauthorized" pulling a public image**: run `docker logout` first, then pull
  without credentials.

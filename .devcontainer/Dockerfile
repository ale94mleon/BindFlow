FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# Not pinned: Ubuntu drops old package versions from its archive, so pins would break the build.
# hadolint ignore=DL3008
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       git wget curl ca-certificates bzip2 sudo build-essential pandoc openssh-client \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN useradd -m -s /bin/bash user \
    && echo "user ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/user \
    && chmod 0440 /etc/sudoers.d/user

# micromamba: the frontend used by install/installer.sh (static binary, no base env)
RUN ARCH="$(uname -m)" \
    && case "$ARCH" in x86_64) MM_PLAT=linux-64 ;; aarch64) MM_PLAT=linux-aarch64 ;; *) echo "Unsupported arch $ARCH" && exit 1 ;; esac \
    && curl -fsSL "https://micro.mamba.pm/api/micromamba/${MM_PLAT}/latest" | tar -xj -C /usr/local bin/micromamba

USER user
WORKDIR /home/user

ENV MAMBA_ROOT_PREFIX=/home/user/micromamba \
    MAMBA_EXE=/usr/local/bin/micromamba

# The BindFlow and BindFlow-TOFF environments are created by install/installer.sh in
# postCreateCommand (.devcontainer/post-create.sh), because the editable install needs the
# mounted repository. The BindFlow env goes first on PATH so it also works without
# activation (VS Code tasks, debugger); BindFlow finds BindFlow-TOFF next to it.
ENV PATH="/home/user/micromamba/envs/BindFlow/bin:/home/user/.local/bin:${PATH}"

# Interactive shells: same env, properly activated (runs its activate.d scripts)
RUN micromamba shell init -s bash -r "$MAMBA_ROOT_PREFIX" \
    && echo "micromamba activate BindFlow 2>/dev/null || true" >> ~/.bashrc

# Claude Code (native installer, installs into /home/user/.local/bin)
RUN curl -fsSL https://claude.ai/install.sh | bash

# Created here so the login volume mounted on it is owned by 'user'
RUN mkdir -p /home/user/.claude

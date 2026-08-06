# GITHUB_PUSH.md —— 项目推送交接说明

> 给新窗口 / 新会话的 AI 助手：按本文件操作即可直接提交并推送，无需重新配置。
> 用法：在对话框粘贴本文件内容（或说一句"按 GITHUB_PUSH.md 执行提交推送"）。

## 项目与仓库信息

| 项 | 值 |
|---|---|
| 项目路径 | `/home/jack/gpcr-ptm` |
| 远程仓库 | https://github.com/okjackishere/gpcr-ptm |
| 远程地址 | `git@github.com:okjackishere/gpcr-ptm.git`（SSH，已配置） |
| 默认分支 | `main`（已跟踪 `origin/main`） |
| git 身份 | `okjackishere <okjackishere@users.noreply.github.com>`（本地已配置） |
| SSH 密钥 | `~/.ssh/id_ed25519`（公钥已登记在 GitHub，免密推送） |

以上均一次性配置完成，持久保存在本机，新会话无需重复设置。

## 提交并推送（每次修改后）

```bash
cd /home/jack/gpcr-ptm
git add -A
git commit -m "对本次改动的简短说明"
git push
```

- 有改动就三条；无改动时 `git push` 会提示 `Everything up-to-date`，属正常。
- `git push` 只上传"已 commit"的内容；改了文件必须先用 `git add` + `git commit`。

## 不要提交的内容（`.gitignore` 已自动排除）

`venv/`、`__pycache__/`、`*.pyc`、生成的 `*_ptm.json` / `*_report.html` / `*_verification.md`、`data/pmid_abstract_cache.json` 等。用 `git status --short` 确认暂存清单里没有这些再提交。

## 故障排查

```bash
ssh -T git@github.com        # 期望输出: Hi okjackishere! You've successfully authenticated...
git remote -v                # 期望: origin  git@github.com:okjackishere/gpcr-ptm.git (fetch/push)
git status -sb               # 期望: ## main...origin/main
git config user.name         # 期望: okjackishere
git config user.email        # 期望: okjackishere@users.noreply.github.com
```

若身份丢失，重新设置（仅对本仓库）：
```bash
cd /home/jack/gpcr-ptm
git config user.name "okjackishere"
git config user.email "okjackishere@users.noreply.github.com"
```

## 换新电脑（而非新窗口）时

新电脑没有这把 SSH 密钥，二选一：
- **复制密钥**：把 `~/.ssh/id_ed25519` 与 `id_ed25519.pub` 拷到新电脑，然后 `git clone git@github.com:okjackishere/gpcr-ptm.git`
- **或 HTTPS**：`git clone https://github.com/okjackishere/gpcr-ptm.git`，推送时用 GitHub Personal Access Token 认证

## 给 AI 的一句话交代

> 项目在 `/home/jack/gpcr-ptm`，改完代码后请 `git add -A && git commit -m "说明" && git push` 推送到 GitHub（SSH 已配置好，可直接推送）。

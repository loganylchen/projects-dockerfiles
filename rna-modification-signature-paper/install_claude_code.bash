curl -fsSL https://fnm.vercel.app/install | bash

source ~/.bashrc
# 新开一个terminal，让 fnm 生效
fnm install 24.3.0
fnm default 24.3.0
fnm use 24.3.0
 
# 安装 claude-code
npm install -g @anthropic-ai/claude-code --registry=https://registry.npmmirror.com
 
# 初始化配置
node --eval "
    const homeDir = os.homedir(); 
    const filePath = path.join(homeDir, '.claude.json');
    if (fs.existsSync(filePath)) {
        const content = JSON.parse(fs.readFileSync(filePath, 'utf-8'));
        fs.writeFileSync(filePath,JSON.stringify({ ...content, hasCompletedOnboarding: true }, 2), 'utf-8');
    } else {
        fs.writeFileSync(filePath,JSON.stringify({ hasCompletedOnboarding: true }), null, 'utf-8');
    }"

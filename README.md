# lie_algebra

## 構成

- doc
    + lie_algebra.md : 三次元回転と並進のリー代数の説明
- lie_se2.py : 2次元回転と並進
- lie_so3.py : 3次元回転
- lie_se3.py : 3次元回転と並進
- lie_ply.py : 3次元回転と並進 (ウサギの点群の例)

## 実行方法

### lie_ply.py

lie_se3.py と同じ階層に置き、data ディレクトリに bun000.ply を置く。
plyファイルの読み込みには https://github.com/dranjan/python-plyfile を使う。

bun000.ply は https://graphics.stanford.edu/data/3Dscanrep/ より拾ってくる。

### lie_ply.py 以外

そのまま python スクリプトを実行。
適当にいじると点群位置が変わる。

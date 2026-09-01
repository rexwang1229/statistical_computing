let 今天 = new Date();
今天.setHours(0,0,0,0);

let 結束 = new Date(2026,11,23);

let 相差天數 = (結束 - 今天) / (1000 * 60 * 60 * 24);

document.getElementById("倒數").textContent = "距離完成這學期的課程還有" + 相差天數 + "天";


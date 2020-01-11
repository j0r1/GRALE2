#include "cppqt.h"
#include <QPainter>
#include <QUuid>
#include <QVariant>
#include <QTextDocument>
#include <QDebug>
#include <QGraphicsScene>
#include <QKeyEvent>
#include <assert.h>

QString GenerateUuid()
{
    return QUuid::createUuid().toString();
}

PointInfo::~PointInfo()
{
}

TriangleInfo::TriangleInfo(const QString &triangUuid, PointInfo *p0, PointInfo *p1, PointInfo *p2)
	: m_triangleUuid(triangUuid)
{
	m_p[0] = p0;
	m_p[1] = p1;
	m_p[2] = p2;

	for (int i = 0 ; i < 3 ; i++)
	{
		assert(m_p[i]);
		assert(m_p[i]->m_triangles.find(this) == m_p[i]->m_triangles.end());
		m_p[i]->m_triangles.insert(this);
	}
}

TriangleInfo::~TriangleInfo()
{
	for (int i = 0 ; i < 3 ; i++)
	{
		assert(m_p[i]);

		auto it = m_p[i]->m_triangles.find(this);
		assert(it != m_p[i]->m_triangles.end());

		m_p[i]->m_triangles.erase(it);
	}
}

Layer::Layer(const QString &name) 
    : m_name(name), m_transform(1, 0, 0, 1, 0, 0), m_uuid(GenerateUuid())
{
}

Layer::~Layer()
{
	clearAllPoints();
}

QString Layer::addPoint(double x, double y, QString label, double timedelay)
{
	QString uuid = GenerateUuid();
	m_points[uuid] = new PointInfo(uuid, x, y, label, timedelay);
	return uuid;
}

QString Layer::addPoint(QVariant xy, QVariant label, QVariant timedelay)
{
	QString uuid = GenerateUuid();
	setPoint(uuid, xy, label, timedelay);
	return uuid;
}

PointInfo *Layer::setPoint(QString uuid, double x, double y, QString label, double timedelay)
{
	auto it = m_points.find(uuid);
	PointInfo *pInf = nullptr;

	if (it == m_points.end())
	{
		pInf = new PointInfo(uuid, x, y, label, timedelay);
		m_points[uuid] = pInf;
	}
	else
	{
		pInf = it.value();
		assert(pInf);

		pInf->m_x = x;
		pInf->m_y = y;
		pInf->m_label = label;
		pInf->m_timedelay = timedelay;
	}
	return pInf;
}

void Layer::setPoint(QString uuid, QVariant xy, QVariant label, QVariant timedelay)
{
	double x = 0;
	double y = 0;
	QList<QVariant> l = xy.toList();
	if (l.size() >= 2)
	{
		x = l[0].toDouble();
		y = l[1].toDouble();
	}

	QString labelStr;
	if (!label.isNull())
		labelStr = label.toString();

	double timedelayDouble = qQNaN();
	if (!timedelay.isNull())
		timedelayDouble = timedelay.toDouble();

	setPoint(uuid, x, y, labelStr, timedelayDouble);
}

QVariant Layer::clearPoint(QString uuid)
{
	QMap<QString, QVariant> affectedTriangles;

	auto pIt = m_points.find(uuid);
	if (pIt == m_points.end())
		return affectedTriangles;

	PointInfo *pInf = pIt.value();
	assert(pInf);

	// These triangles will be removed, remove them from other points

	// Destructor will also remove iterators from this list!
	auto triangles = pInf->m_triangles;
	for (auto tInf : triangles)
	{
		QString triangKey = tInf->m_triangleUuid;
		QList<QVariant> triangPoints { tInf->m_p[0]->m_uuid, tInf->m_p[1]->m_uuid, tInf->m_p[2]->m_uuid };
		affectedTriangles[triangKey] = triangPoints;

		auto tIt = m_triangles.find(triangKey);
		assert(tIt != m_triangles.end());
		m_triangles.erase(tIt);

		delete tInf; // deletes itself from the lists of the points
	}
	delete pInf;
	m_points.erase(pIt);

	return affectedTriangles;
}

void Layer::clearAllPoints()
{
	for (auto t : m_triangles)
		delete t;
	for (auto p : m_points)
		delete p;

    m_points.clear();
    m_triangles.clear();
}

PointInfo *Layer::getPointInfo(QString uuid)
{
	auto it = m_points.find(uuid);
	if (it == m_points.end())
		return nullptr;
	return it.value();
}

inline QVariant pointInfoToVariant(double x, double y, const QString &label, double timedelay)
{
	QMap<QString, QVariant> values;
	QList<QVariant> xy;

	xy.push_back(x);
	xy.push_back(y);
	values["xy"] = xy;

	if (label.isNull())
		values["label"] = QVariant();
	else
		values["label"] = label;
	if (!qIsNaN(timedelay))
		values["timedelay"] = timedelay;
	else
		values["timedelay"] = QVariant();
	return values;
}

inline QVariant pointInfoToVariant(const PointInfo &ptInf)
{
	return pointInfoToVariant(ptInf.m_x, ptInf.m_y, ptInf.m_label, ptInf.m_timedelay);
}

QVariant Layer::getPoint(QString uuid) const
{
	auto it = m_points.find(uuid);
	if (it == m_points.end())
		return QVariant();

	assert(it.value());
	return pointInfoToVariant(*it.value());
}

QVariant Layer::getPoints(bool transformed)
{
	QMap<QString, QVariant> points;
	for (auto it = m_points.begin() ; it != m_points.end() ; ++it)
	{
		if (!transformed)
			points[it.key()] = pointInfoToVariant(*it.value());
		else
		{
			auto x = it.value()->m_x;
			auto y = it.value()->m_y;

			auto xyTransformed = m_transform.map(QPointF(x, y));
			points[it.key()] = pointInfoToVariant(xyTransformed.x(), xyTransformed.y(), it.value()->m_label, it.value()->m_timedelay);
		}
	}
	return points;
}

QVariant Layer::getTrianglesContainingPoint(QString uuid)
{
    auto it = m_points.find(uuid);
    if (it == m_points.end())
        return QList<QVariant>();

	QList<QVariant> x;
	for (auto i : it.value()->m_triangles)
		x.push_back(i->m_triangleUuid);

    return x;
}

QString Layer::addTriangle(const QString &p0Key, const QString &p1Key, const QString &p2Key, const QString &triangUuid)
{
	TriangleInfo *pTriang = addTriangleInfo(p0Key, p1Key, p2Key, triangUuid);
	return pTriang->m_triangleUuid;
}

TriangleInfo *Layer::addTriangleInfo(const QString &p0Key, const QString &p1Key, const QString &p2Key, const QString &triangUuid)
{
    QString uuid = triangUuid;

    if (uuid.isEmpty())
        uuid = GenerateUuid();
	else
	{
		// Check that triangle doesn't exist
		assert(m_triangles.find(uuid) == m_triangles.end());
	}

	assert(p0Key != p1Key);
	assert(p1Key != p2Key);
	assert(p0Key != p2Key);
	PointInfo *p0Inf = getPointInfo(p0Key);
	PointInfo *p1Inf = getPointInfo(p1Key);
	PointInfo *p2Inf = getPointInfo(p2Key);
	assert(p0Inf);
	assert(p1Inf);
	assert(p2Inf);
    
	TriangleInfo *pTriangInf = new TriangleInfo(uuid, p0Inf, p1Inf, p2Inf);
	m_triangles[uuid] = pTriangInf;

    return pTriangInf;
}

QVariant Layer::getTriangles() const
{
	QMap<QString, QVariant> triangs;

	for (auto it = m_triangles.begin() ; it != m_triangles.end() ; ++it)
	{
		const TriangleInfo &tInf = *(it.value());
		QList<QVariant> pts { tInf.m_p[0]->m_uuid, tInf.m_p[1]->m_uuid, tInf.m_p[2]->m_uuid };
		triangs[it.key()] = pts;
	}
    return triangs;
}

const TriangleInfo *Layer::getTriangleInfo(const QString &triangUuid) const
{
	auto it = m_triangles.find(triangUuid);
	if (it == m_triangles.end())
		return nullptr;
	return it.value();
}

QVariant Layer::getTriangle(const QString &triangUuid) const
{
	auto it = m_triangles.find(triangUuid);
	if (it == m_triangles.end())
		return QVariant();

	const TriangleInfo &tInf = *(it.value());
	QList<QVariant> pts { tInf.m_p[0]->m_uuid, tInf.m_p[1]->m_uuid, tInf.m_p[2]->m_uuid };
	return pts;
}

void Layer::clearTriangle(const QString &triangUuid)
{
	auto it = m_triangles.find(triangUuid);
	if (it == m_triangles.end())
		return;

	const TriangleInfo *tInf = it.value();
	delete tInf;

	m_triangles.erase(it);
}

EmptyGraphicsItem::EmptyGraphicsItem(bool childrenBoundingRect, QGraphicsItem *pParent) 
	: m_childrenBoundingRect(childrenBoundingRect), QGraphicsItem(pParent)
{
}

EmptyGraphicsItem::~EmptyGraphicsItem()
{
}

QRectF EmptyGraphicsItem::boundingRect() const
{
	if (!m_childrenBoundingRect)
		return QRectF();

	auto r = QRectF(0,0,0,0);
	for (auto &i : childItems())
	{
		if (i->isVisible())
			r |= i->mapToParent(i->boundingRect()).boundingRect();
	}
	return r;
}

void EmptyGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
}

QBrush TriangleItem::s_normalBrush(QColor(192, 192, 192, 128));
QBrush TriangleItem::s_selectedBrush(QColor(192, 192, 255, 128));
QPen TriangleItem::s_normalPen(QBrush(Qt::darkGray), 2);
QPen TriangleItem::s_selectedPen(QBrush(Qt::darkGray), 2);
bool TriangleItem::s_brushInit = false;

// Just to have some default
TriangleItem::TriangleItem(QGraphicsItem *pParent) 
    : QGraphicsItem(pParent), m_pLayer(nullptr)
{
}

TriangleItem::TriangleItem(const QString &uuid, Layer *pLayer, const QMap<QString, QPointF> &pointCoords, QGraphicsItem *pParent) 
    : QGraphicsItem(pParent),
      m_triangleUuid(uuid), m_pLayer(pLayer), m_ptsCoords(pointCoords)
{
    if (!s_brushInit)
    {
        s_brushInit = true;
        s_normalPen.setCosmetic(true);
        s_selectedPen.setCosmetic(true);
    }
    
    updatePolygon();
    setFlag(QGraphicsItem::ItemIsSelectable);
}

TriangleItem::~TriangleItem()
{
}

bool TriangleItem::contains(const QPointF &p) const
{
    return m_polygon.containsPoint(p, Qt::OddEvenFill);
}

void TriangleItem::updatePointPosition(const QString &ptUuid, QPointF &pos)
{
    auto it = m_ptsCoords.find(ptUuid);
    assert(it != m_ptsCoords.end());
    it.value() = pos;

    updatePolygon();
}

void TriangleItem::updatePolygon()
{
    QVector<QPointF> pts;
    for (auto p : m_ptsCoords)
        pts.push_back(p);
    assert(pts.size() == 3);

	m_polygon = QPolygonF(pts);
	prepareGeometryChange();
}

void TriangleItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
	if (isSelected())
	{
		painter->setPen(s_selectedPen);
		painter->setBrush(s_selectedBrush);
	}
	else
	{
		painter->setPen(s_normalPen);
		painter->setBrush(s_normalBrush);
	}
	painter->drawPolygon(m_polygon);
}

QRectF TriangleItem::boundingRect() const
{
	// TODO: handle the cosmetic pen better!
	QRectF r = m_polygon.boundingRect();
	r.adjust(-0.2, -0.2, 0.2, 0.2);
	return r;
}

TriangleItem *TriangleItem::getTriangleItem(QGraphicsItem *pItem)
{
	return dynamic_cast<TriangleItem *>(pItem);
}

SyncTextEditItem::SyncTextEditItem(Layer *pLayer, const QString &uuid, const QString &label, QGraphicsItem *pParent)
    : m_pLayer(pLayer), m_pointUuid(uuid), QGraphicsTextItem(label, pParent)
{
}

SyncTextEditItem::~SyncTextEditItem()
{
}

void SyncTextEditItem::keyReleaseEvent(QKeyEvent *pEvt)
{
    QGraphicsTextItem::keyReleaseEvent(pEvt);
    syncValueAndCenter();

	if (pEvt->key() == Qt::Key_Escape)
		clearFocus();
}

void SyncTextEditItem::inputMethodEvent(QInputMethodEvent *pEvt)
{
    QGraphicsTextItem::inputMethodEvent(pEvt);
    syncValueAndCenter();
}

void SyncTextEditItem::focusOutEvent(QFocusEvent *pEvt)
{
    QGraphicsTextItem::focusOutEvent(pEvt);
    syncValueAndCenter();

    const PointInfo *pPtInf = m_pLayer->getPointInfo(m_pointUuid);
    assert(pPtInf);

    const QString &label = pPtInf->m_label;
    if (label.isEmpty())
        setVisible(false);
}

void SyncTextEditItem::recenter()
{
    auto r = boundingRect();
    auto dx = (r.right()+r.left())/2.0;
    setTransform(QTransform(0.1,0,0,-0.1,0,0));
    setTransform(QTransform(1,0,0,1,-dx, 10), true);
}

void SyncTextEditItem::syncValueAndCenter()
{
    QString label = document()->toPlainText();
    const PointInfo *pPtInf = m_pLayer->getPointInfo(m_pointUuid);
    assert(pPtInf);

    PointInfo ptOld(*pPtInf);
    m_pLayer->setPoint(m_pointUuid, ptOld.m_x, ptOld.m_y, label, ptOld.m_timedelay);

    SceneBase *pScene = dynamic_cast<SceneBase *>(scene());
    if (pScene)
        pScene->onPointLabelChanged(m_pLayer->getUuid(), m_pointUuid, ptOld.m_label, label);

    recenter();
}

LayerObjectGraphicsItem::LayerObjectGraphicsItem()
	: m_pLayer(nullptr), EmptyGraphicsItem(true, nullptr)
{
}

LayerObjectGraphicsItem::LayerObjectGraphicsItem(Layer *pLayer, const QString &objectUuid, QGraphicsItem *pParent)
    : m_pLayer(pLayer), m_uuid(objectUuid), EmptyGraphicsItem(true, pParent)
{
}

LayerObjectGraphicsItem::~LayerObjectGraphicsItem()
{
}

QBrush PointGraphicsItemBase::s_fontBrush(Qt::green);

PointGraphicsItemBase::PointGraphicsItemBase()
	: LayerObjectGraphicsItem(nullptr, "", nullptr)
{
	m_pTxt = nullptr;
	m_pTdTxt = nullptr;
}

PointGraphicsItemBase::PointGraphicsItemBase(Layer *pLayer, const QString &uuid, const QString &label, double timedelay, QGraphicsItem *pParent)
    : LayerObjectGraphicsItem(pLayer, uuid, pParent)
{
    setMovable();

    m_pTxt = new SyncTextEditItem(pLayer, uuid, label, this);
    m_pTxt->setFlag(QGraphicsItem::ItemIsFocusable);
    m_pTxt->setDefaultTextColor(QColor(Qt::cyan));
    m_pTxt->setTextInteractionFlags(Qt::TextEditorInteraction);
    m_pTxt->recenter();

    if (label.isEmpty())
        m_pTxt->setVisible(false);

    m_pTdTxt = new QGraphicsSimpleTextItem(this);
    m_pTdTxt->setPen(QPen(s_fontBrush, 0.2));
    m_pTdTxt->setBrush(s_fontBrush);

    if (!qIsNaN(timedelay))
    {
        setTDText(timedelay);
        m_pTdTxt->setVisible(true);
    }
    else
        m_pTdTxt->setVisible(false);

    SceneBase *pScene = dynamic_cast<SceneBase *>(scene());
    if (pScene)
        setTransform(pScene->getPointTransform());
}

PointGraphicsItemBase::~PointGraphicsItemBase()
{
}

PointGraphicsItemBase *PointGraphicsItemBase::getPointGraphicsItem(QGraphicsItem *pItem)
{
	PointGraphicsItemBase *pPointItem = dynamic_cast<PointGraphicsItemBase *>(pItem);
	return pPointItem;
}

bool PointGraphicsItemBase::isNormalPoint() const
{
	return (dynamic_cast<const SinglePointGraphicsItem *>(this) != nullptr);
}

bool PointGraphicsItemBase::isMovable() const
{
    auto f = flags();
    if (f&QGraphicsItem::ItemIsMovable && f&QGraphicsItem::ItemIsSelectable)
        return true;
    return false;
}

void PointGraphicsItemBase::setMovable(bool v)
{
    setFlag(QGraphicsItem::ItemIsSelectable, v);
    setFlag(QGraphicsItem::ItemIsMovable, v);
}

void PointGraphicsItemBase::fetchSettings()
{
    Layer *pLayer = getLayer();
    const PointInfo *pPtInf = pLayer->getPointInfo(getUuid());
    assert(pPtInf);

	fetchSettings(*pPtInf);
}

void PointGraphicsItemBase::fetchSettings(const PointInfo &pPtInf)
{
    Layer *pLayer = getLayer();
    QPointF pos = pLayer->getImageTransform().map(QPointF(pPtInf.m_x, pPtInf.m_y));
    setPos(pos);

    if (pPtInf.m_label.isEmpty())
        m_pTxt->setVisible(false);
    else
    {
        m_pTxt->setPlainText(pPtInf.m_label);
        m_pTxt->setVisible(true);
        m_pTxt->recenter();
    }

    if (pPtInf.isTimedelayValid())
    {
        m_pTdTxt->setVisible(true);
        setTDText(pPtInf.m_timedelay);
    }
    else
        m_pTdTxt->setVisible(false);
}
        
void PointGraphicsItemBase::setTDText(double timedelay)
{
	QString s = QString::number(timedelay) + " days";
	m_pTdTxt->setText(s);

	auto r = m_pTdTxt->boundingRect();
	auto dx = (r.right()+r.left())/2.0;
	m_pTdTxt->setTransform(QTransform(0.1,0,0,-0.1,0,0));
	m_pTdTxt->setTransform(QTransform(1,0,0,1,-dx, -20), true);
}

QVariant PointGraphicsItemBase::syncPosition()
{
    Layer *pLayer = getLayer();
	QPointF pos0 = pos();
    auto p = pLayer->getImageTransform().inverted().map(pos0);

    const PointInfo *pPtInf = pLayer->getPointInfo(getUuid());
    assert(pPtInf);

	// Map to scene coords
	QPointF oldXY = pLayer->getImageTransform().map(QPointF(pPtInf->m_x, pPtInf->m_y));

    QList<QVariant> oldPos { oldXY.x(), oldXY.y() };
    QList<QVariant> newPos { pos0.x(), pos0.y() }; // also screen coords
    QList<QVariant> oldAndNewPos { oldPos, newPos };

	// Save the pixel coords
    pLayer->setPoint(getUuid(), p.x(), p.y(), pPtInf->m_label, pPtInf->m_timedelay);

    return oldAndNewPos;
}

QVariant PointGraphicsItemBase::itemChange(GraphicsItemChange change, const QVariant &value)
{
    if (change == QGraphicsItem::ItemSelectedChange)
        onSelected(value.toBool());
    else if (change == QGraphicsItem::ItemPositionHasChanged)
    {
        LayerGraphicsItemBase *pLayerItem = static_cast<LayerGraphicsItemBase *>(parentItem());
        assert(pLayerItem);

        pLayerItem->onPointPositionChanged(getUuid(), value.toPointF());
    }

    return LayerObjectGraphicsItem::itemChange(change, value);
}

void PointGraphicsItemBase::toggleFocus()
{
    if (m_pTxt->hasFocus())
        m_pTxt->clearFocus();
    else
    {
        m_pTxt->setFocus();
        m_pTxt->setVisible(true);
    }
}

QBrush SinglePointGraphicsItem::s_normalBrush(QColor(255, 255, 0, 128));
QBrush SinglePointGraphicsItem::s_selectedBrush(QColor(0, 0, 255, 128));
QPen SinglePointGraphicsItem::s_normalPen(QBrush(Qt::gray), 0.2);
QPen SinglePointGraphicsItem::s_selectedPen(QBrush(Qt::darkGray), 0.2);
QPen SinglePointGraphicsItem::s_pointSelNormalPen(QBrush(Qt::red), 0.5);
QPen SinglePointGraphicsItem::s_pointSelSelectedPen(QBrush(Qt::blue), 0.5);
bool SinglePointGraphicsItem::s_init = false;

SinglePointGraphicsItem::SinglePointGraphicsItem(Layer *pLayer, const QString &uuid, const QString &label, double timedelay, bool isPointSelect, QGraphicsItem *pParent)
    : PointGraphicsItemBase(pLayer, uuid, label, timedelay, pParent), 
	  m_isPointSelect(isPointSelect)
{
	if (!s_init)
	{
		s_init = true;
		s_pointSelSelectedPen.setCosmetic(true);
		s_pointSelNormalPen.setCosmetic(true);
	}
}

SinglePointGraphicsItem::~SinglePointGraphicsItem()
{
}

QRectF SinglePointGraphicsItem::boundingRect() const
{
	QRectF r = PointGraphicsItemBase::boundingRect();
	r |= QRectF(-1.2,-1.2,2.4,2.4);
	return r;
}

void SinglePointGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
	if (isSelected())
	{
		auto &pen = (m_isPointSelect)?s_pointSelSelectedPen:s_selectedPen;
		painter->setPen(pen);
		painter->setBrush(s_selectedBrush);
	}
	else
	{
		auto &pen = (m_isPointSelect)?s_pointSelNormalPen:s_normalPen;
		painter->setPen(pen);
		painter->setBrush(s_normalBrush);
	}

	if (!m_isPointSelect)
	{
		painter->drawEllipse(QRectF(-1, -1, 2, 2));

		painter->setPen(s_normalPen);
		painter->drawLine(QPointF(-0.2,-0.2), QPointF(0.2, 0.2));
		painter->drawLine(QPointF(-0.2, 0.2), QPointF(0.2, -0.2));
	}
	else
	{
		painter->drawLine(QPointF(-1,-1), QPointF(1, 1));
		painter->drawLine(QPointF(-1, 1), QPointF(1, -1));
	}
}

void SinglePointGraphicsItem::onSelected(bool v)
{
	/*
    if (v)
    {
        m_pCircle->setPen(s_selectedPen);
        m_pCircle->setBrush(s_selectedBrush);
    }
    else
    {
        m_pCircle->setPen(s_normalPen);
        m_pCircle->setBrush(s_normalBrush);
    }*/
}

QPen MatchPointGraphicsItemBase::s_normalLinePen(QBrush(Qt::red), 0.2);
QPen MatchPointGraphicsItemBase::s_selectedLinePen(QBrush(Qt::blue), 0.2);
QBrush MatchPointGraphicsItemBase::s_rectBrush(QColor(0, 0, 0, 1));
QPen MatchPointGraphicsItemBase::s_rectPen(MatchPointGraphicsItemBase::s_rectBrush, 0.01);

MatchPointGraphicsItemBase::MatchPointGraphicsItemBase(Layer *pLayer, const QString &uuid, const QString &label, QGraphicsItem *pParent)
	: PointGraphicsItemBase(pLayer, uuid, label, qQNaN(), pParent)
{
	m_pRect = new QGraphicsRectItem(-1.1, -1.1, 2.2, 2.2, this);
	m_pRect->setPen(s_rectPen);
	m_pRect->setBrush(s_rectBrush);
}

MatchPointGraphicsItemBase::~MatchPointGraphicsItemBase()
{
}

MatchPointGraphicsItemCross::MatchPointGraphicsItemCross(Layer *pLayer, const QString &uuid, const QString &label, QGraphicsItem *pParent)
	: MatchPointGraphicsItemBase(pLayer, uuid, label, pParent)
{
	m_pLine1 = new QGraphicsLineItem(-1,-1, 1, 1, this);
	m_pLine2 = new QGraphicsLineItem(-1, 1, 1,-1, this);
	onSelected(false);
}

MatchPointGraphicsItemCross::~MatchPointGraphicsItemCross()
{
}

void MatchPointGraphicsItemCross::onSelected(bool v)
{
	if (v)
	{
		m_pLine1->setPen(s_selectedLinePen);
		m_pLine2->setPen(s_selectedLinePen);
	}
	else
	{
		m_pLine1->setPen(s_normalLinePen);
		m_pLine2->setPen(s_normalLinePen);
	}
}

MatchPointGraphicsItemCircle::MatchPointGraphicsItemCircle(Layer *pLayer, const QString &uuid, const QString &label, QGraphicsItem *pParent)
	: MatchPointGraphicsItemBase(pLayer, uuid, label, pParent)
{
	m_pCircle = new QGraphicsEllipseItem(-1, -1, 2, 2, this);
	onSelected(false);
}

MatchPointGraphicsItemCircle::~MatchPointGraphicsItemCircle()
{
}

void MatchPointGraphicsItemCircle::onSelected(bool v)
{
	if (v)
		m_pCircle->setPen(s_selectedLinePen);
	else
		m_pCircle->setPen(s_normalLinePen);
}

void SceneBase::setPointTransform(const QTransform &t)
{
	m_pointTransform = t;

	for (auto item : items())
	{
		PointGraphicsItemBase *pItem = dynamic_cast<PointGraphicsItemBase *>(item);
		if (pItem)
			pItem->setTransform(t);
	}
}

LayerGraphicsItemBase::LayerGraphicsItemBase(Layer *pLayer, PointType pType, bool childrenBoundingRect, QGraphicsItem *pParent)
	: m_pLayer(pLayer), m_pointType(pType), EmptyGraphicsItem(childrenBoundingRect, pParent)
{
	m_pointsVisible = true;
}

LayerGraphicsItemBase::~LayerGraphicsItemBase()
{
}

void LayerGraphicsItemBase::onPointPositionChanged(const QString &pointUuid, QPointF pos)
{
	const PointInfo *pInf = m_pLayer->getPointInfo(pointUuid);
	if (!pInf)
		return;

	const QSet<TriangleInfo *> &triangles = pInf->m_triangles;
	for (auto t : triangles)
	{
		QString tUuid = t->m_triangleUuid;
		assert(m_triangleItems.find(tUuid) != m_triangleItems.end());
		m_triangleItems[tUuid]->updatePointPosition(pointUuid, pos);
	}
}

PointGraphicsItemBase *LayerGraphicsItemBase::createPoint(const PointInfo &pt)
{
	PointGraphicsItemBase *pItem = nullptr;
	switch(m_pointType)
	{
	case Normal:
		pItem = new SinglePointGraphicsItem(m_pLayer, pt.m_uuid, pt.m_label, pt.m_timedelay, false, this);
		break;
	case PointSelect:
		pItem = new SinglePointGraphicsItem(m_pLayer, pt.m_uuid, pt.m_label, pt.m_timedelay, true, this);
		break;
	case Cross:
		pItem = new MatchPointGraphicsItemCross(m_pLayer, pt.m_uuid, pt.m_label, this);
		break;
	case Circle:
	default:
		pItem = new MatchPointGraphicsItemCircle(m_pLayer, pt.m_uuid, pt.m_label, this);
	}
	assert(m_pLayer);

	QPointF pos = m_pLayer->getImageTransform().map(QPointF(pt.m_x, pt.m_y));
	pItem->setZValue(2);

	// WARNING! This triggers position change event, need to make sure that
	//          the onPointPositionChanged callback doesn't cause a loop this
	//          way
	pItem->setPos(pos);
	pItem->setVisible(m_pointsVisible);

	return pItem;
}

TriangleItem *LayerGraphicsItemBase::createTriangle(const TriangleInfo *triangInf)
{
	QMap<QString, QPointF> pointCoords;

	assert(m_pLayer);
	assert(triangInf);
	assert(triangInf->m_p[0]);
	assert(triangInf->m_p[1]);
	assert(triangInf->m_p[2]);

	for (int i = 0 ; i < 3 ; i++)
	{
		const PointInfo &p = *(triangInf->m_p[i]);
		pointCoords[triangInf->m_p[i]->m_uuid] = m_pLayer->getImageTransform().map(QPointF(p.m_x, p.m_y));
	}

	TriangleItem *pItem = new TriangleItem(triangInf->m_triangleUuid, m_pLayer, pointCoords, this);
	pItem->setZValue(1);

	assert(m_triangleItems.find(triangInf->m_triangleUuid) == m_triangleItems.end());
	return pItem;
}

template<class T>
inline void deleteUnmarkedItems(QMap<QString, T *> &items)
{
	for (auto it = items.begin() ; it != items.end() ; /* */ )
	{
		T *pItem = it.value();
		if (!pItem->isMarked())
		{
			delete pItem;
			it = items.erase(it);
		}
		else
			++it;
	}
}

void LayerGraphicsItemBase::updatePoints()
{
	for (auto p : m_pointItems)
		p->clearMarker();
	for (auto t : m_triangleItems)
		t->clearMarker();

	QTransform trans(1, 0, 0, 1, 0, 0);
	SceneBase *pScene = dynamic_cast<SceneBase *>(scene());
	if (pScene)
		trans = pScene->getPointTransform();

	auto points = m_pLayer->getPointInfos();
	for (auto it = points.begin(); it != points.end() ; ++it)
	{
		QString uuid = it.key();
		PointInfo &pt = *(it.value());

		PointGraphicsItemBase *pItem = nullptr;
		auto it2 = m_pointItems.find(uuid);

		if (it2 == m_pointItems.end())
		{
			pItem = createPoint(pt);
			m_pointItems[uuid] = pItem;
		}
		else
		{
			pItem = it2.value();
			pItem->fetchSettings(pt);
			pItem->setVisible(m_pointsVisible);
		}
		pItem->setTransform(trans);
		pItem->mark();
	}

	auto triangles = m_pLayer->getTriangleInfos();
	for (auto it = triangles.begin() ; it != triangles.end() ; ++it)
	{
		const QString &triangUuid = it.key();
		TriangleInfo *triangInf = it.value();

		TriangleItem *pItem = nullptr;
		auto it2 = m_triangleItems.find(triangUuid);
		if (it2 == m_triangleItems.end())
		{
			pItem = createTriangle(triangInf);
			m_triangleItems[triangUuid] = pItem;
		}
		else
			pItem = it2.value();

		pItem->mark();
	}

	deleteUnmarkedItems<TriangleItem>(m_triangleItems);
	deleteUnmarkedItems<PointGraphicsItemBase>(m_pointItems);
}

void LayerGraphicsItemBase::setPointsVisible(bool v) 
{ 
	m_pointsVisible = v;
	updatePoints();
}

void LayerGraphicsItemBase::updatePointList(const QList<QVariant> &points, const QList<QVariant> &triangles)
{
	// TODO
	updatePoints();
}

void LayerGraphicsItemBase::setLabel(const QString &uuid, const QString &label)
{
	PointInfo *pInf = m_pLayer->getPointInfo(uuid);
	assert(pInf);
	pInf->m_label = label;

	auto it = m_pointItems.find(uuid);
	assert(it != m_pointItems.end());
	PointGraphicsItemBase *pItem = it.value();
	assert(pItem);

	pItem->fetchSettings(*pInf);
}

void LayerGraphicsItemBase::movePoint(const QString &uuid, double x, double y)
{
	PointInfo *pInf = m_pLayer->getPointInfo(uuid);
	assert(pInf);
	
	QPointF transPos = m_pLayer->getImageTransform().inverted().map(QPointF(x, y));
	pInf->m_x = transPos.x();
	pInf->m_y = transPos.y();
	
	auto it = m_pointItems.find(uuid);
	assert(it != m_pointItems.end());
	PointGraphicsItemBase *pItem = it.value();
	assert(pItem);

	pItem->setPos(x, y);

	// Also update the triangles
	const QSet<TriangleInfo *> &triangles = pInf->m_triangles;
	auto pointUuid = pInf->m_uuid;
	QPointF p(x, y);

	for (auto t : triangles)
	{
		QString tUuid = t->m_triangleUuid;
		assert(m_triangleItems.find(tUuid) != m_triangleItems.end());

		m_triangleItems[tUuid]->updatePointPosition(pointUuid, p);
	}
}

QString LayerGraphicsItemBase::addPoint(double x, double y, const QString &label, double timedelay, const QString &uuid)
{
	QString ptUuid = uuid;

	if (uuid.isEmpty())
		ptUuid = GenerateUuid();
	else
	{
		assert(m_pLayer->getPointInfo(uuid) == nullptr);
	}

	QPointF transPos = m_pLayer->getImageTransform().inverted().map(QPointF(x, y));
	PointInfo *pInf = m_pLayer->setPoint(ptUuid, transPos.x(), transPos.y(), label, timedelay);

	assert(m_pointItems.find(ptUuid) == m_pointItems.end());
	PointGraphicsItemBase *pItem = createPoint(*pInf);
	m_pointItems[ptUuid] = pItem;

	QTransform trans(1, 0, 0, 1, 0, 0);
	SceneBase *pScene = dynamic_cast<SceneBase *>(scene());
	if (pScene)
		trans = pScene->getPointTransform();
	pItem->setTransform(trans);

	return ptUuid;
}

void LayerGraphicsItemBase::setPoint(const QString &uuid, double x, double y, const QString &label, double timedelay)
{
	PointInfo *pInf = m_pLayer->getPointInfo(uuid);
	assert(pInf);

	pInf->m_x = x;
	pInf->m_y = y;
	pInf->m_label = label;
	pInf->m_timedelay = timedelay;

	PointGraphicsItemBase *pItem = nullptr;
	auto it = m_pointItems.find(uuid);
	if (it == m_pointItems.end())
	{
		pItem = createPoint(*pInf);
		m_pointItems[uuid] = pItem;
	}
	else
	{
		pItem = it.value();
		pItem->fetchSettings(*pInf);

		// Also update the triangles
		const QSet<TriangleInfo *> &triangles = pInf->m_triangles;
		auto pointUuid = pInf->m_uuid;
		QPointF p(x, y);

		for (auto t : triangles)
		{
			QString tUuid = t->m_triangleUuid;
			assert(m_triangleItems.find(tUuid) != m_triangleItems.end());

			m_triangleItems[tUuid]->updatePointPosition(pointUuid, p);
		}
	}
	QTransform trans(1, 0, 0, 1, 0, 0);
	SceneBase *pScene = dynamic_cast<SceneBase *>(scene());
	if (pScene)
		trans = pScene->getPointTransform();
	pItem->setTransform(trans);
}

QString LayerGraphicsItemBase::addTriangle(const QString &p0, const QString &p1, const QString &p2, const QString &triangleUuid)
{
	QString tUuid = triangleUuid;

	if (triangleUuid.isEmpty())
		tUuid = GenerateUuid();

	assert(m_pLayer->getTriangleInfo(tUuid) == nullptr);

	TriangleInfo *pInf = m_pLayer->addTriangleInfo(p0, p1, p2, tUuid);

	auto it = m_triangleItems.find(tUuid);
	assert(it == m_triangleItems.end());
	m_triangleItems[tUuid] = createTriangle(pInf);
	return tUuid;
}

void LayerGraphicsItemBase::clearTriangle(const QString &triangleUuid)
{
	m_pLayer->clearTriangle(triangleUuid);

	auto it = m_triangleItems.find(triangleUuid);
	if (it != m_triangleItems.end())
	{
		delete it.value(); // delete the item
		m_triangleItems.erase(it); // remove the entry
	}
}

QVariant LayerGraphicsItemBase::clearPoint(const QString &uuid)
{
	PointInfo *pInf = m_pLayer->getPointInfo(uuid);
	if (!pInf)
		return QVariant();

	auto pIt = m_pointItems.find(uuid);
	if (pIt != m_pointItems.end())
	{
		delete pIt.value();
		m_pointItems.erase(pIt);
	}

	// Remove the point from the layer data
	QVariant affectedTriangles = m_pLayer->clearPoint(uuid);
	
	auto triangList = affectedTriangles.toMap();

	// Remove the triangle items
	for (auto i = triangList.begin() ; i != triangList.end() ; ++i)
	{
		QString triangUuid = i.key();

		auto it = m_triangleItems.find(triangUuid);
		if (it != m_triangleItems.end())
		{
			delete it.value(); // delete the item
			m_triangleItems.erase(it); // remove the entry
		}
	}

	return affectedTriangles;
}

PointGraphicsItemBase *LayerGraphicsItemBase::getPointItem(const QString &uuid)
{
	auto pIt = m_pointItems.find(uuid);
	if (pIt == m_pointItems.end())
		return nullptr;

	return pIt.value();
}

QVariant splitPointsAndTriangles(Layer &layer, bool ignoreRemainingPoints)
{
	auto points = layer.getPointInfos();
	auto triangles = layer.getTriangleInfos();

	QList<QVariant> imageInfo;

	QSet<PointInfo *> allPoints;
	QSet<TriangleInfo *> allTriangles;

	for (auto p : points)
		allPoints.insert(p);
	for (auto t : triangles)
		allTriangles.insert(t);

	while (true)
	{
		if (allTriangles.size() == 0)
			break;

		QSet<TriangleInfo *> curTriangles;
		QSet<PointInfo *> curPoints;
	
		TriangleInfo *pStartTriang = *allTriangles.begin();
		curTriangles.insert(pStartTriang);

		auto tIt = allTriangles.find(pStartTriang);
		if (tIt == allTriangles.end())
			return QString("Internal error: couldn't find triangle in allTriangles");
		allTriangles.erase(tIt);

		for (int i = 0 ; i < 3 ; i++)
		{
			PointInfo *pPointInf = pStartTriang->m_p[i];
			curPoints.insert(pPointInf);

			auto pIt = allPoints.find(pPointInf);
			if (pIt == allPoints.end())
				return QString("Internal error: couldn't find point in allPoints");
			allPoints.erase(pIt);
		}

		/*
		qDebug() << "Initial points:" << endl;

		for (auto p : curPoints)
			qDebug() << "    " << p->m_uuid;

		qDebug() << "";
		*/
	
		// Look for other triangles in these points:
		while (true)	
		{
			bool foundNew = false;

			QSet<PointInfo *> newPoints;

			for (auto p : curPoints)
			{
				//qDebug() << "Checking triangles in " << p->m_uuid;
				for (auto k : p->m_triangles)
				{
					// Check that we didn't process this triangle yet
					auto tIt = allTriangles.find(k);
					if (tIt != allTriangles.end())
					{

						foundNew = true;
						curTriangles.insert(k);
						allTriangles.erase(tIt);

						for (int i = 0 ; i < 3 ; i++)
							newPoints.insert(k->m_p[i]);
					}
				}
			}

			if (!foundNew)
				break;

			curPoints.unite(newPoints);
			allPoints.subtract(newPoints);
		}

		/*
		qDebug() << "Final points:" << endl;
		for (auto p : curPoints)
			qDebug() << "    " << p->m_uuid;

		qDebug() << "";
		*/

		QMap<QString, QVariant> inf;

		QMap<QString, QVariant> pointsVar;
		QMap<QString, QVariant> triangsVar;
		for (auto p : curPoints)
			pointsVar[p->m_uuid] = pointInfoToVariant(p->m_x, p->m_y, p->m_label, p->m_timedelay);

		for (auto t : curTriangles)
		{
			QList<QVariant> pts { t->m_p[0]->m_uuid, t->m_p[1]->m_uuid, t->m_p[2]->m_uuid };
			triangsVar[t->m_triangleUuid] = pts;
		}

		inf["points"] = pointsVar;
		inf["triangles"] = triangsVar;

		imageInfo.push_back(inf);
	}

	if (!ignoreRemainingPoints && !allPoints.empty())
	{
		QList<QVariant> pointsLeft;

		for (auto p : allPoints)
			pointsLeft.push_back(p->m_uuid);

		QMap<QString, QVariant> err;

		err["error"] = QString("Can't split points layer in multiple images: still points left after considering all triangulations");
		err["pointsleft"] = pointsLeft;
		return err;
	}

    return imageInfo;
}

